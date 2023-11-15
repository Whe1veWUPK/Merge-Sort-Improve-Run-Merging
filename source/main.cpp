#include<chrono>
#include<iostream>
#include<thread>
#include<vector>
#include<queue>
#include<deque>
#include<fstream>
#include<cmath>
#include "LoserTree.hpp"
#include "FileOperator.hpp"
#include "Timer.hpp"
//int curRunNum
int currentRunNum = 0;
//vector that store the each run's length
std::vector<int>runLength;
//for get the best mergesequence
struct Node {
    int parent;
    int leftChild;
    int rightChild;
    int weight;
};
//read data from the disk
void readFromDisk(Buffer* inputBuffer, FileOperator* fileOperator, int startPosition, std::string filePath) {
    //get inputBuffer's Size 
    int predictSize = inputBuffer->getBufferSize();
    //
    inputBuffer->resize(predictSize);
    //read info from file
    fileOperator->writeToInputBuffer(filePath, predictSize, startPosition, inputBuffer);
}
//generate Run
void generateRun(LoserTree* loserTree, Buffer* inputBuffer, Buffer* outputBuffer,int*overLen,bool*isReset) {
    *overLen = -1;
    *isReset = false;
    //get two buffers' size
    int bufferSize = inputBuffer->getCurSize();
    //num of valid data
    int count = 0;
    for (int i = 0; i < bufferSize; ++i) {
       //get data from tree
        outputBuffer->append(i, loserTree->getTop().getValue());
        // move data into the tree
        loserTree->adjust(loserTree->getTop().getKey(), inputBuffer->buffer[i]);
        ++count;
        //if all the data in the tree is invalid
        if (loserTree->needToCreateNewRun()) {
            //std::cout << loserTree->getTop().getValue() << " \n";
            loserTree->setValid();
            //record the overLen
            *overLen = count;
            *isReset = true;
        }
        
    }
    //if the tree has not set the data's validity
    if (*overLen == -1) {
        //count is equal to the bufferSize
        *overLen = count;
    }
}
//write data to the disk
void writeToDisk(Buffer* outputBuffer, FileOperator* FileOperator, bool isFirst, std::string  filePath) {

    //get outputBuffer's size
    int predictSize = outputBuffer->getCurSize();
    if (outputBuffer->Empty()) {
        return;
    }
    //write the buffer's data to disk
    FileOperator->writeToFile(filePath, predictSize, isFirst, outputBuffer);
}
//read the run's data that will exhaust first
void readRunData(Buffer*freeBuffer,FileOperator*fileOperator,std::string filePath,std::deque<int>*workingBufferQueue,int numOfMergeWays,int runNum,int*runCurPosition,int*runStartPosition,bool isFirst,bool*isOver,int*curRunNum,int startPos,int dataSize,std::queue<int>*freeIndex,int bufferSize){
    //the process has been done in intial part
    if(isFirst){
        return;
    }
    char runPath[128];
    //find the run that will exhaust first
    int key = -1;
    int minValue = (((1 << 30) - 1) * 2 + 1);
    for (int i = 0; i < numOfMergeWays;++i){
        if(isOver[i]){
            // if the run is over, then just skip
            continue;
        }
        int curSize = freeBuffer[workingBufferQueue[i].back()].getCurSize();
        if(freeBuffer[workingBufferQueue[i].back()].buffer[curSize-1]<minValue){
            key = i;
            minValue = freeBuffer[workingBufferQueue[i].back()].buffer[curSize - 1];
        }
    }
    //if all runs have been read
    if(key==-1){
        return;
    }
    //write data to the buffer
    //get the new Buffer's index
    int newKey = freeIndex->front();
    freeIndex->pop();
    //add the new Key to the working buffer queue
    workingBufferQueue[key].push_back(newKey);
    //update run's cur pos
    if(runCurPosition[startPos+key]+bufferSize<=runLength[startPos+key]){
        sprintf(runPath, "Run\\Run%d.txt", startPos + key);
        fileOperator->writeToInputBuffer(runPath, bufferSize, runCurPosition[startPos+key], &freeBuffer[newKey]);
        runCurPosition[startPos+key] += bufferSize;

    }
    else if((runCurPosition[startPos+key]+bufferSize>runLength[startPos+key])&&(runCurPosition[startPos+key]<runLength[startPos+key])){
        sprintf(runPath, "Run\\Run%d.txt", startPos + key);
        fileOperator->writeToInputBuffer(runPath, runLength[startPos+key] - runCurPosition[startPos+key],runCurPosition[startPos+key],&freeBuffer[newKey]);
        runCurPosition[startPos+key] = runLength[startPos+key];

    }
    //set the over info and curRunNum info
    if(runCurPosition[startPos+key]>=runLength[startPos+key]){
        //the run is over
        if(isOver[key]){
            return;
        }
        else{
            isOver[key] = true;
            --*curRunNum;
        }
    }
    //std::cout << freeBuffer[workingBufferQueue[key].back()].getCurSize() << "\n";
    int curSize=freeBuffer[workingBufferQueue[key].back()].getCurSize();
    if(curSize==0){
       //the run is over
       //1.resize the buffer
        freeBuffer[workingBufferQueue[key].back()].resize(bufferSize);
       //2. pop_back it
        workingBufferQueue[key].pop_back();
       //3  push it to the freeIndex
        freeIndex->push(newKey);
    }
}
//k-way-merge
void startMerge(LoserTree*loserTree,Buffer*outputBuffer,Buffer*freeBuffer,std::deque<int>*workingBufferQueue,int*runNumNotMerged,bool*isOver,std::queue<int>*freeIndex){
    // get bufferSize
    int bufferSize = outputBuffer->getBufferSize();
    //resize outputbuffer
    outputBuffer->resize(bufferSize);
    if(*runNumNotMerged==0){
        //no run,then over 
        return;
    }
    // Do until the outputBuffer becomes full or all runs are over
    for (int i = 0; i < bufferSize;++i){
        if(*runNumNotMerged==0){
            //no run, then return
            return;
        }
        //write the winner into the output buffer
        outputBuffer->append(i, loserTree->getTop().getValue());
        //get the winner's key and adjust the loserTree
        int key = loserTree->getTop().getKey();
        if(workingBufferQueue[key].size()==0&&isOver[key]){
            //if the run is Over
            std::cout << "The queue is empty and the run is over\n";
            loserTree->adjust(-1, false);
            continue;
        }
        else{
            //the run is not Over
            // record the index
            int indexTemp = workingBufferQueue[key].front();
            //store the buffer
            Buffer temp(bufferSize);
            for (int i = 0; i <freeBuffer[indexTemp].getCurSize() ;++i){
                temp.append(i, freeBuffer[indexTemp].buffer[i]);
            }
            //set location 
            temp.setCurLoation(freeBuffer[indexTemp].getCurLocation());

            //std::cout << "The buffer's size is:" << freeBuffer[workingBufferQueue[key].front()].getCurSize() << "\n";


            //move curPos
            freeBuffer[workingBufferQueue[key].front()].moveCurPos();
            if(freeBuffer[workingBufferQueue[key].front()].isOver()){
                //the buffer's data has all been merged
                //1.free the buffer
                freeBuffer[workingBufferQueue[key].front()].resize(bufferSize);

                //2.move to the next buffer in the queue
                workingBufferQueue[key].pop_front();
                if(workingBufferQueue[key].size()==0&&isOver[key]){
                    // the run is over
                    //add invalid value to the loserTree
                    loserTree->adjust(-1,false);
                    std::cout << "The Run " << key << " is Over\n";
                    //change the num
                    --*runNumNotMerged;
                    continue;
                }
                else if(workingBufferQueue[key].size()==0&&!isOver[key]){
                    //the run is not Over
                    //exit the merge, to wait the I/O, recovery the data
                    //1. back the output buffer
                    outputBuffer->backCurPos();
                    //2. recovery the data
                    freeBuffer[indexTemp] = Buffer(bufferSize);
                    for (int i = 0; i < temp.getCurSize();++i){
                        freeBuffer[indexTemp].append(i, temp.buffer[i]);
                    }
                    freeBuffer[indexTemp].setCurLoation(temp.getCurLocation());
                    //std::cout << "the copy buffer " << freeBuffer[indexTemp].getCurSize() << "\n";
                    // 3. recover the queue
                    workingBufferQueue[key].push_back(indexTemp);
                    //4. exit
                    return;
                }
                //3. add new value to the loserTree
                int curLoc = freeBuffer[workingBufferQueue[key].front()].getCurLocation();
                int newValue = freeBuffer[workingBufferQueue[key].front()].buffer[curLoc];
                loserTree->adjust(key, newValue);
                //4. update freeIndex
                freeIndex->push(indexTemp);
            }
            else{
                //the buffer is not over
                //add new value to the loserTree
                int curLoc = freeBuffer[workingBufferQueue[key].front()].getCurLocation();
                int newValue = freeBuffer[workingBufferQueue[key].front()].buffer[curLoc];
                loserTree->adjust(key, newValue);
            }
        }

    }
}
void createRuns(std::string filePath){
    FileOperator file;
    std::fstream ifs;
    ifs.open(filePath, std::ios::in);
    std::string line;
    char dataOutputPath[128];
    int curPos = 0;
    int curNum = 0;
    std::fstream ofs;
    sprintf(dataOutputPath, "Run\\Run%d.txt", curNum);
    ofs.open(dataOutputPath,std::ios::out);
    while(std::getline(ifs,line)){
        if(curPos<runLength[curNum]){
            ofs << line;
            ofs << "\n";
            ++curPos;

        }
        else{
            ofs.close();
            curPos = 0;
            ++curNum;
            sprintf(dataOutputPath, "Run\\Run%d.txt", curNum);
            ofs.open(dataOutputPath,std::ios::out);
            ofs << line;
            ofs << "\n";
            ++curPos;
        }
    }
    ifs.close();
}
void startSort(int bufferSize,int treeSize) {
    FileOperator file;
    int dataNum = file.getDataSize("Input.txt");
    //output the data size
    std::cout << "The data size is " << dataNum << "\n";

    Buffer* buffers = new Buffer[3];
    for (int i = 0; i < 3; ++i) {
        buffers[i] = Buffer(bufferSize);
    }

    FileOperator  fileOperator;
    std::string filepath = "Input.txt";
    std::string outputPath = "Output.txt";
    int dataSize = fileOperator.getDataSize(filepath);
    int curLocation = 0;
    //initial loserTree
    Buffer tempBuffer = Buffer(treeSize);
    fileOperator.writeToInputBuffer(filepath, treeSize, curLocation, &tempBuffer);
    LoserTree* loserTree = new LoserTree(tempBuffer.buffer, treeSize);
    //update curLocation
    curLocation += treeSize;
    bool isFirst = true;
    int pin = 0;
    //initial one buffer
    fileOperator.writeToInputBuffer(filepath,buffers[0].getBufferSize(), curLocation, &buffers[0]);
    //update curLocation
    curLocation += buffers[0].getCurSize();
    int curLen = 0;
    // tag to refer the cur Run whether is over
    int overLen = 0;
    //loop until the whole file has been read
    //std::cout << "initial buffer already\n";
    while (curLocation < dataSize) {
        //std::cout << "CurLocation is: " << curLocation << "\n";
        //read from the disk
        std::thread t1(readFromDisk,&buffers[(pin+1)%3],&fileOperator,curLocation,filepath);
        //run generation
        overLen = 0;
        bool isReset=false;
        std::thread t2(generateRun, loserTree, &buffers[pin], &buffers[(pin+2)%3],&overLen,&isReset);
        t2.join();
        //std::cout << "Run generated successfully\n";
        //write the output buffer to the disk
        std::thread t3(writeToDisk, &buffers[(pin + 2) % 3], &fileOperator, isFirst, outputPath);
        // if the run is over
        if (overLen<buffers[(pin+2)%3].getCurSize()||isReset) {
            //std::cout << overLen << "<"<< buffers[(pin + 2) % 3].getCurSize()<<"\n";

            //compute the curLen
            curLen += overLen;
            //record the cur run's lenth
            //std::cout << "CurLen is:" << curLen << "\n";
            runLength.emplace_back(curLen);
            //reset the curLen
            curLen = buffers[(pin + 2) % 3].getCurSize() - overLen;
            
        }
        else {
            curLen += buffers[(pin + 2) % 3].getCurSize();
        }
        t1.join();
       //update cur location
        curLocation += buffers[(pin + 1) % 3].getCurSize();
        
        t3.join();
        //update First flag
        if (isFirst) {
            isFirst = false;
        }
        //clear buffer
        buffers[pin].resize(buffers[pin].getBufferSize());
        buffers[(pin + 2) % 3].resize(buffers[pin].getBufferSize());
       //update next pin
        pin = (pin + 1) % 3;

    }
    //there is a buffer that has already read data from disk but has not been write to loserTree
    overLen = 0;
    bool isReset = false;
    generateRun(loserTree, &buffers[pin], &buffers[(pin + 2)%3], &overLen,&isReset);
    writeToDisk(&buffers[(pin + 2) % 3], &fileOperator, false, outputPath);
    if (overLen < buffers[(pin + 2) % 3].getCurSize()) {
        //std::cout << overLen << "<" << buffers[(pin + 2) % 3].getCurSize() << "\n";

        //compute the curLen
        curLen += overLen;
        //record the cur run's lenth
        //std::cout << "CurLen is:" << curLen << "\n";
        runLength.emplace_back(curLen);
        //reset the curLen
        curLen = buffers[(pin + 2) % 3].getCurSize() - overLen;

    }
    else {
        curLen += buffers[(pin + 2) % 3].getCurSize();
    }
    //std::cout <<"The loser Tree's size is: " <<loserTree->getSize() << "\n";
    //the loser tree has remaining data
    //input  maxValue  into the tree
    int maxValue = (1<< 30);
    int maxNum = 0;
    Buffer finalBuffer(loserTree->getSize());

    while (maxNum<loserTree->getSize()) {
        finalBuffer.append(maxNum, loserTree->getTop().getValue());
        loserTree->adjust (maxValue,false);
        ++curLen;
        //if all data in the loser tree is invalid
        if (loserTree->needToCreateNewRun()) {
            
            loserTree->setValid();
            //record the run's length
            runLength.emplace_back(curLen);
            
            curLen = 0;
        }
        ++maxNum;
    }
    if(curLen!=0){
       //record the final run's length
       runLength.emplace_back(curLen);
       //write the final buffer to the disk

    }
    writeToDisk(&finalBuffer, &fileOperator,false,"Output.txt");
    createRuns("Output.txt");
}
//get Best Merge Sequence
void getBestMergeSequence() {
    int cost = 0;
    std::cout << "Best merge Sequence are as followed\n";
    //get runLength's size
    int size = runLength.size();
    std::vector<Node> nodes;
    //initial nodes
    for (int i = 0; i < size; ++i) {
        Node node;
        node.weight = runLength[i];
        //-1 means that it has only one node
        node.parent = -1;
        node.leftChild = -1;
        node.rightChild = -1;
        //append the node
        nodes.emplace_back(node);
    }
    //std::cout << "Runs' length are as followed: \n";

    //for (int i = 0; i < size; ++i) {
    //    std::cout << runLength[i] << " ";
    //}
    //std::cout << " \n";
    for (int i = 0; i < size-1; ++i) {
        for (int i = 0; i < nodes.size(); ++i) {
            if (nodes[i].parent == -1) {
                std::cout << nodes[i].weight << " ";
            }

        }
        std::cout << "\n";
        int min1 = (1 << 30);
        int index1 = 0;
        int min2 = (1 << 30);
        int index2 = 0;
        //find two min value
        for (int i = 0; i < nodes.size(); ++i) {
            if (nodes[i].parent == -1 && nodes[i].weight < min1) {
                min1 = nodes[i].weight;
                index1 = i;
            }
        }
        for (int i = 0; i < nodes.size(); ++i) {
            if (nodes[i].parent == -1 && nodes[i].weight < min2 && i != index1) {
                min2 = nodes[i].weight;
                index2 = i;
            }
        }
        //merge and create new node
        std::cout << "Merge " << nodes[index1].weight << " and " << nodes[index2].weight << "\n";
        Node newNode;
        newNode.weight = nodes[index1].weight + nodes[index2].weight;
        cost += newNode.weight;
        newNode.parent = -1;
        newNode.leftChild = index1;
        newNode.rightChild = index2;
        nodes.emplace_back(newNode);
        nodes[index1].parent = nodes.size() - 1;
        nodes[index2].parent = nodes.size() - 1;
    }
    std::cout << "The merge cost is: " << cost << "\n";
}
void mergeRuns(int k,int bufferSize){
    //the input and output file path
    std::string inputFilePath = "Input.txt";
    std::string outputFilePath = "Output.txt";
    char runPath[128];
    //the vector that store the file path
    std::vector<std::string> filePath;
    filePath.push_back(inputFilePath);
    filePath.push_back(outputFilePath);
    FileOperator f;
    //get the total data size
    int dataSize = 0;
    for (int i = 0; i < runLength.size();++i){
        dataSize += runLength[i];
    }
    int runNum = runLength.size();
    int runNotFinishedMergeNum = runNum; //record the num of runs that have not been merged
    std::cout << "The data size is: " << dataSize << "\n";
    std::cout << "The number of run is: " << runNum << "\n";

    //initiate buffers and the loserTree
    Buffer *freeBuffer = new Buffer[2 * k];
    //cur pin
    int curPin = 0;
    for (int i = 0; i < 2 * k;++i){
        freeBuffer[i] = Buffer(bufferSize);
       
    }
    Buffer *outputBuffer = new Buffer[2];
    for (int i = 0; i < 2;++i){
        outputBuffer[i] = Buffer(bufferSize);
    }
    std::cout << "Output Buffer's Size is: \n";
    for (int i = 0; i < 2;++i){
        std::cout << outputBuffer[i].getBufferSize() << " ";
    }
    std::cout << "\n";
    //intiate the free queue
    std::queue<int>freeIndex;

    for (int i = 0; i < 2 * k;++i){
        freeIndex.push(i);
    }

    std::cout << freeIndex.front() << "\n";
    //initiate the working buffer's queue
    std::deque<int>*workingBufferQueue=new std::deque<int>[k];
    // get each run's initial start Position
    int *runStartPosition = new int[runNum];
    for (int i = 0; i < runNum;++i){
       if(i==0){
        runStartPosition[i] = 0;
       }
       else{
        runStartPosition[i] = runStartPosition[i - 1] + runLength[i - 1];
       }
    }
    //initial each run's curPosition
    int *runCurPosition = new int[runNum];
    for (int i = 0; i < runNum; ++i){
        runCurPosition[i] = 0;
    }
    std::cout << "runStartPosition: \n";
    for (int i = 0; i < runNum;++i){
       std::cout << runStartPosition[i] << " ";
    }
    std::cout << "\n";
    std::cout << "runCurPosition: \n";
    for (int i = 0; i < runNum;++i){
       std::cout << runCurPosition[i] << " ";
    }
    std::cout << "\n";
    if (runNum <= k){
      //if num of run is less than or equal to k
      //flag to report the run is over
       bool *isOver = new bool[runNum];
       for (int i = 0; i < runNum;++i){
        isOver[i] = false;
       }
       // initiate working buffer
       for (int i = 0; i < runNum;++i){
        int index = freeIndex.front();
        workingBufferQueue[i].push_back(index);
        freeIndex.pop();
        
       }
       //read data from disk to working buffer
        for (int i = 0; i < runNum; ++i){
          if (runCurPosition[i] + bufferSize <= runLength[i]){
              sprintf(runPath, "Run\\Run%d.txt", i);
              f.writeToInputBuffer(runPath, bufferSize, runCurPosition[i], &freeBuffer[workingBufferQueue[i].front()]);
              runCurPosition[i] += bufferSize;
            }
          else if (runCurPosition[i] + bufferSize > runLength[i]){
              sprintf(runPath, "Run\\Run%d.txt", i);
              f.writeToInputBuffer(runPath, runLength[i] - runCurPosition[i], runCurPosition[i], &freeBuffer[workingBufferQueue[i].front()]);
              runCurPosition[i] = runLength[i];
            }
        }
       //initiate data
       int *data = new int[runNum];
       for (int i = 0; i < runNum;++i){
          data[i] = freeBuffer[workingBufferQueue[i].front()].buffer[0];
       }

    //    for (int i = 0; i < runNum;++i){
    //       std::cout << data[i] << " ";
    //    }
    //    std::cout << "\n";
       //build the loserTree
       LoserTree loserTree(data, runNum);
       //find the run that will exhaust first
       int key = -1;
       int minValue = (((1 << 30)-1)*2+1);
       //std::cout << minValue << "\n";
       for (int i = 0; i < runNum;++i){
          int curSize = freeBuffer[workingBufferQueue[i].back()].getCurSize();
          if (freeBuffer[workingBufferQueue[i].back()].buffer[curSize-1]<minValue){
              minValue = freeBuffer[workingBufferQueue[i].back()].buffer[curSize - 1];
              key = i;
          }
       }
       int index = freeIndex.front();
       freeIndex.pop();
       workingBufferQueue[key].push_back(index);
       //write data to the buffer

       if(runCurPosition[key]+bufferSize<=runLength[key]){
           sprintf(runPath, "Run\\Run%d.txt", key);
           f.writeToInputBuffer(runPath, bufferSize, runCurPosition[key], &freeBuffer[index]);
           runCurPosition[key] += bufferSize;            
       }
       else if(runCurPosition[key]+bufferSize>runLength[key]){
            sprintf(runPath, "Run\\Run%d.txt", key);
            f.writeToInputBuffer(runPath, runLength[key] - runCurPosition[key],runCurPosition[key],&freeBuffer[index]);
            runCurPosition[key] = runLength[key];
       }
       int curSize = freeBuffer[workingBufferQueue[key].back()].getCurSize();
       if(curSize==0){
                   //the run is over
                   //1.resize the buffer
                   freeBuffer[workingBufferQueue[key].back()].resize(bufferSize);
                   //2. pop_back it
                   workingBufferQueue[key].pop_back();
                   //3  push it to the freeIndex
                   freeIndex.push(index);
        }
       //set activeOutputBuffer
       int activeOutputBuffer = 0;
       bool isFirst = true;
    //    for (int i = 0; i < runNum;++i){
    //         std::cout << "workingBuffer " << (i) << "'s size is: " << workingBufferQueue[i].size() << "\n";
    //    }
       //startMerge
       std::fstream fout;
       if (isFirst) {
        //如果是某一轮次的第一个 outputBuffer输出 则重写文档
        fout.open(outputFilePath, std::ios::out | std::ios::trunc);
       }
       else {
        
        fout.open(outputFilePath, std::ios::out | std::ios::app);
       }
       fout.close();
       int totalRunNum = runNum;
       while((runNotFinishedMergeNum>0)){
        for (int i = 0; i < totalRunNum;++i){
              std::cout << "Run " << i << "'s data: \n";
              for (int j = 0; j <workingBufferQueue[i].size();++j){
                    std::cout << "FreeBuffer " << workingBufferQueue[i].at(j) <<"\n";
                    for (int k = 0; k < freeBuffer[workingBufferQueue[i].at(j)].getCurSize();++k){
                        std::cout << freeBuffer[workingBufferQueue[i].at(j)].buffer[k] << " ";
                    }
                    std::cout << "\n";
              }
              std::cout << "\n";
        }
        
        // thread 1: do the k-way-merge in memory
        std::thread t1(startMerge, &loserTree, &outputBuffer[activeOutputBuffer], freeBuffer, workingBufferQueue,&runNotFinishedMergeNum,isOver, &freeIndex);
        // thread 2: write the another output buffer's data to the disk
        std::thread t2(writeToDisk, &outputBuffer[(activeOutputBuffer+1)%2], &f, isFirst, outputFilePath);
        // wait thread 1 to finish
        t1.join();
        // thread 3: read the run's data that it will first exhaust
        std::thread t3(readRunData, freeBuffer, &f, outputFilePath, workingBufferQueue, totalRunNum, totalRunNum, runCurPosition, runStartPosition, isFirst,isOver,&runNum,0,dataSize,&freeIndex, bufferSize);
        // wait thread 2 and thread 3 to finish
        t2.join();
        t3.join();
        isFirst = false;
        // std::cout << "Output buffer generated!\n";
        // test = false;
        activeOutputBuffer = (activeOutputBuffer + 1) % 2;
        }
       // write the last output buffer into the disk
        writeToDisk(&outputBuffer[(activeOutputBuffer + 1) % 2], &f, false, outputFilePath);
    }
    else{
        //calculate the merge pass round
        double passRound = (std::ceil(std::log((long double)runNum) / std::log((long double)k)));
        std::cout << "Need " << passRound << " rounds\n";
        int curRoundNum = runNum;
        int curRun = 0;
        // set the curRound Merge's file
        int writeIndex = 0;
        int readIndex = 1;
        //runPath
        char runPath[128];
        std::vector<int> nextRoundRunLength;
        //do the Merge
        while((passRound--)>0){

          //the number that next round
          int nextRoundNum = 0;
          //clear the nextRoundRunLength vector
          nextRoundRunLength.clear();
          bool isFirst = true;
          //initial txt
          std::fstream fout;
          if (isFirst) {
           //如果是某一轮次的第一个 outputBuffer输出 则重写文档
            fout.open(outputFilePath, std::ios::out | std::ios::trunc);
          }
          else {
           fout.open(outputFilePath, std::ios::out | std::ios::app);
          }
          fout.close();
          //std::cerr << "K is: " << k << "\n";
          // start the round Merge
          int curRoundTotalNum = curRoundNum;
          while(curRoundNum>0){
          //std::cerr << "K is: " << k << "\n"; 
                //resize the freeBuffer
                for (int i = 0; i < 2 * k;++i){
                      freeBuffer[i].resize(bufferSize);
                }
                //resize the outputBuffer
                for (int i = 0; i < 2;++i){
                      outputBuffer[i].resize(bufferSize);
                }
                //clear the freeIndex
                while(!freeIndex.empty()){
                      freeIndex.pop();
                }
                for (int i = 0; i < 2 * k;++i){
                      freeIndex.push(i);
                }
              if(curRoundNum>=k){         
                int runLen = 0;
                // k - way - merge
                //  flag to record if the run is read over
                bool *isOver = new bool[k];
                for (int i = 0; i < k; ++i){
                        isOver[i] = false;
                }
                // initiate working buffer
                for (int i = 0; i < k;++i){
                    int index = freeIndex.front();
                    workingBufferQueue[i].push_back(index);
                    freeIndex.pop();
                }
                //read data from disk to working buffer
                for (int i = 0; i < k; ++i){
                     if (runCurPosition[curRun+i] + bufferSize <= runLength[curRun+i]){
                         sprintf(runPath, "Run\\Run%d.txt", curRun+i);
                         f.writeToInputBuffer(runPath, bufferSize, runCurPosition[curRun +i], &freeBuffer[workingBufferQueue[i].front()]);
                         runCurPosition[curRun +i] += bufferSize;
                     }
                     else if (runCurPosition[curRun+i]+bufferSize>runLength[curRun+i]&&(runCurPosition[curRun+i]<runLength[curRun+i])){
                         sprintf(runPath, "Run\\Run%d.txt", curRun + i);
                         f.writeToInputBuffer(runPath, runLength[curRun+i]-runCurPosition[curRun+i], runCurPosition[curRun +i], &freeBuffer[workingBufferQueue[i].front()]);
                         runCurPosition[curRun + i] = runLength[curRun + i];
                        }

                }
                //initiate data
                //std::cerr << "K is: " << k << "\n";
                int*treeData = new int[k];
                //std::cout << "K is: " << k << "\n";
                for (int i = 0; i < k;++i){
                     treeData[i] = freeBuffer[workingBufferQueue[i].front()].buffer[0];
                }
                std::cout << "Initial data\n";
                for (int i = 0; i < k;++i){
                     std::cout << treeData[i] << " ";
                }
                std::cout << "\n";
                // build the loserTree
                LoserTree loserTree(treeData, k);
                std::cout << "Loser Tree Winner: " << loserTree.getTop().getValue() << "\n";
                //find the first run that will exhaust first
                int key = -1;
                int minValue = (((1 << 30)-1)*2+1);
                //std::cout << minValue << "\n";
                for (int i = 0; i < k;++i){
                int curSize = freeBuffer[workingBufferQueue[i].back()].getCurSize();
                   if (freeBuffer[workingBufferQueue[i].back()].buffer[curSize-1]<minValue){
                       minValue = freeBuffer[workingBufferQueue[i].back()].buffer[curSize - 1];
                       key = i;
                    }
                }
                int index = freeIndex.front();
                freeIndex.pop();
                workingBufferQueue[key].push_back(index);
                //write data to the buffer

                if(runCurPosition[curRun+key]+bufferSize<=runLength[curRun+key]){
                    sprintf(runPath, "Run\\Run%d.txt", key + curRun);
                    f.writeToInputBuffer(runPath, bufferSize, runCurPosition[curRun + key], &freeBuffer[index]);
                    runCurPosition[curRun + key] += bufferSize;            
                }
                else if(runCurPosition[curRun+key]+bufferSize>runLength[curRun+key]&&(runCurPosition[curRun+key]<runLength[curRun+key])){
                     sprintf(runPath, "Run\\Run%d.txt", key + curRun);
                     f.writeToInputBuffer(runPath,runLength[curRun+key] - runCurPosition[curRun+key],runCurPosition[curRun+key],&freeBuffer[index]);
                     runCurPosition[curRun + key] = runLength[curRun + key];
                }
                int curSize = freeBuffer[workingBufferQueue[key].back()].getCurSize();
                if(curSize==0){
                   //the run is over
                   //1.resize the buffer
                   freeBuffer[workingBufferQueue[key].back()].resize(bufferSize);
                   //2. pop_back it
                   workingBufferQueue[key].pop_back();
                   //3  push it to the freeIndex
                   freeIndex.push(index);
                }
                //set activeOutputBuffer
                int activeOutputBuffer = 0;
                int totalRunNum = k;
                int runNotMergeFinished = k;
                int curRunNum = k;
                //int testNum = 1111;
                while((runNotMergeFinished>0)){

                    //  for (int i = 0; i < k;++i){
                    //    std::cout << "Run " << i << "\n";
                    //    for (int j = 0; j < workingBufferQueue[i].size();++j){
                    //        std::cout << "Buffer " << workingBufferQueue[i].at(j) << "\n";
                    //        for (int z = 0; z < freeBuffer[workingBufferQueue[i].at(j)].getCurSize();++z){
                    //            std::cout << freeBuffer[workingBufferQueue[i].at(j)].buffer[z] << " ";
                    //        }
                    //        std::cout << "\n";
                    //    }
                    //  }
                    // thread 1: do the k-way-merge in memory
                    std::thread t1(startMerge, &loserTree, &outputBuffer[activeOutputBuffer], freeBuffer, workingBufferQueue, &runNotMergeFinished, isOver, &freeIndex);
                     // thread 2: write the another output buffer's data to the disk
                     std::thread t2(writeToDisk, &outputBuffer[(activeOutputBuffer+1)%2], &f, isFirst, outputFilePath);
                     // wait thread 1 to finish
                     t1.join();
                     //std::cout << "Merge Success\n";
                     // thread 3: read the run's data that it will first exhaust
                     std::thread t3(readRunData, freeBuffer, &f, inputFilePath, workingBufferQueue, k,curRoundTotalNum, runCurPosition, runStartPosition, isFirst, isOver, &curRunNum, curRun, dataSize, &freeIndex, bufferSize);
                     // wait thread 2 and thread 3 to finish
                     t2.join();
                     //std::cout << "Write Success\n";
                     t3.join();
                    // std::cout << "Read Success\n";
                     // update runLen
                     runLen += outputBuffer[(activeOutputBuffer + 1) % 2].getCurSize();
                     isFirst = false;
                     activeOutputBuffer = (activeOutputBuffer + 1) % 2;
                }
                // write the last output buffer into the disk
                writeToDisk(&outputBuffer[(activeOutputBuffer + 1) % 2], &f, false, outputFilePath);
                runLen += outputBuffer[(activeOutputBuffer + 1) % 2].getCurSize();
                //update curRoundNum
                curRoundNum -= k;
                //update curRun
                curRun += k;
                //add the runLen to the vector
                //std::cout << "RunLen is: " <<runLen << "\n";
                nextRoundRunLength.push_back(runLen);
                ++nextRoundNum;
                //reset runLen
                runLen = 0;
                //free the memory
                delete[] treeData;
              }
              else{
                 //curRoundNum - way - merge
                 //record the run's length
                int runLen = 0;
                // flag to record if the run is read over
                bool *isOver = new bool[curRoundNum];
                for (int i = 0; i < curRoundNum; ++i){
                        isOver[i] = false;
                }
                // initiate working buffer
                for (int i = 0; i < curRoundNum;++i){
                    int index = freeIndex.front();
                    workingBufferQueue[i].push_back(index);
                    freeIndex.pop();
                }

                //read data from disk to working buffer
                for (int i = 0; i < curRoundNum; ++i){
                     if (runCurPosition[curRun+i] + bufferSize <= runLength[curRun + i]){
                         sprintf(runPath, "Run\\Run%d.txt", curRun + i);
                         f.writeToInputBuffer(runPath, bufferSize, runCurPosition[curRun +i], &freeBuffer[workingBufferQueue[i].front()]);
                         runCurPosition[curRun +i] += bufferSize;
                     }
                     else if (runCurPosition[curRun+i]+bufferSize>runLength[curRun+i]&&(runCurPosition[curRun+i]<runLength[curRun+i])){
                         sprintf(runPath, "Run\\Run%d.txt", curRun + i);
                         f.writeToInputBuffer(runPath,runLength[curRun+i]-runCurPosition[curRun+i] , runCurPosition[curRun +i], &freeBuffer[workingBufferQueue[i].front()]);
                         runCurPosition[curRun + i] = runLength[curRun + i];
                     }
                }
                //initiate data
                //std::cerr << "K is: " << k << "\n";
                int*treeData = new int[curRoundNum];
                //std::cout << "K is: " << k << "\n";
                for (int i = 0; i < curRoundNum;++i){
                     treeData[i] = freeBuffer[workingBufferQueue[i].front()].buffer[0];
                }
                //build the loserTree
                LoserTree loserTree(treeData, curRoundNum);
                //find the first run that will exhaust first
                int key = -1;
                int minValue = (((1 << 30)-1)*2+1);
                for (int i = 0; i < curRoundNum;++i){
                int curSize = freeBuffer[workingBufferQueue[i].back()].getCurSize();
                   if (freeBuffer[workingBufferQueue[i].back()].buffer[curSize-1]<minValue){
                       minValue = freeBuffer[workingBufferQueue[i].back()].buffer[curSize - 1];
                       key = i;
                    }
                }
                int index = freeIndex.front();
                freeIndex.pop();
                workingBufferQueue[key].push_back(index);
                //write data to the buffer

                if(runCurPosition[curRun+key]+bufferSize<=runLength[curRun+key]){
                    sprintf(runPath, "Run\\Run%d.txt", curRun + key);
                    f.writeToInputBuffer(runPath, bufferSize, runCurPosition[curRun + key], &freeBuffer[index]);
                    runCurPosition[curRun + key] += bufferSize;            
                }
                else if(runCurPosition[curRun+key]+bufferSize>runLength[curRun+key]&&(runCurPosition[curRun+key]<runLength[curRun+key])){
                     sprintf(runPath, "Run\\Run%d.txt", curRun + key);
                     f.writeToInputBuffer(filePath[readIndex], runLength[curRun+key ] - runCurPosition[curRun+key],runCurPosition[curRun+key],&freeBuffer[index]);
                     runCurPosition[curRun + key] = runLength[curRun + key];
                }
                int curSize=freeBuffer[workingBufferQueue[key].back()].getCurSize();
                if(curSize==0){
                //the run is over
                //1.resize the buffer
                freeBuffer[workingBufferQueue[key].back()].resize(bufferSize);
                //2. pop_back it
                workingBufferQueue[key].pop_back();
                //3  push it to the freeIndex
                freeIndex.push(index);
                }
                // set activeOutputBuffer
                int activeOutputBuffer = 0;
                int totalRunNum = curRoundNum;
                int runNotMergeFinished = curRoundNum;
                int curRunNum = curRoundNum;
                     while ((runNotMergeFinished > 0))
                     {
                          for (int i = 0; i < curRoundNum;++i){
                            std::cout << "Run " << i << "\n";
                            for (int j = 0; j < workingBufferQueue[i].size();++j){
                                std::cout << "Buffer " << workingBufferQueue[i].at(j) << "\n";
                                for (int z = 0; z < freeBuffer[workingBufferQueue[i].at(j)].getCurSize();++z){
                                    std::cout << freeBuffer[workingBufferQueue[i].at(j)].buffer[z] << " ";
                                }
                                std::cout << "\n";
                            }
                          }
                         // thread 1: do the k-way-merge in memory
                         std::thread t1(startMerge, &loserTree, &outputBuffer[activeOutputBuffer], freeBuffer, workingBufferQueue, &runNotMergeFinished, isOver, &freeIndex);
                         // thread 2: write the another output buffer's data to the disk
                         std::thread t2(writeToDisk, &outputBuffer[(activeOutputBuffer + 1) % 2], &f, isFirst, outputFilePath);
                         // wait thread 1 to finish
                         t1.join();
                         // std::cout << "Successfully merged\n";
                         //  thread 3: read the run's data that it will first exhaust
                         std::thread t3(readRunData, freeBuffer, &f, filePath[readIndex], workingBufferQueue, totalRunNum, curRoundTotalNum, runCurPosition, runStartPosition, isFirst, isOver, &curRunNum, curRun, dataSize, &freeIndex, bufferSize);
                         // wait thread 2 and thread 3 to finish
                         t2.join();
                         // std::cout << "Successfully write\n";
                         t3.join();
                         // std::cout << "Successfully read\n";
                         //  update runLen
                         runLen += outputBuffer[(activeOutputBuffer + 1) % 2].getCurSize();
                         isFirst = false;
                         activeOutputBuffer = (activeOutputBuffer + 1) % 2;
                }
                // write the last output buffer into the disk
                writeToDisk(&outputBuffer[(activeOutputBuffer + 1) % 2], &f, false, outputFilePath);
                runLen += outputBuffer[(activeOutputBuffer + 1) % 2].getCurSize();
                //update curRun
                curRun += curRoundNum;
                //update curRoundNum
                curRoundNum =0;
                //add the runLen to the vector
                //std::cout << "RunLen is: " <<runLen << "\n";
                nextRoundRunLength.push_back(runLen);
                ++nextRoundNum;
                //reset runLen
                runLen = 0;
                //free the memory
                delete[] treeData;
              }
          }
          if(nextRoundNum==1){
              return;
          }
          //resetCurrentNum
          currentRunNum = 0;
          // update RunLenth
          runLength.clear();
          for (int i = 0; i < nextRoundRunLength.size();++i){
              runLength.push_back(nextRoundRunLength[i]);
          }
          createRuns(outputFilePath);
          // delete the memory
          delete[] runCurPosition;
          //Reset the run's info
          runCurPosition = new int[nextRoundNum];
          for (int i = 0; i < nextRoundNum;++i){
              runCurPosition[i] =0;
          }
          std::cout << "\n";
          //reset the isFirst flag
          isFirst = true;
          //reset the curRun
          curRun = 0;
          //set the curRound num
          curRoundNum = nextRoundNum;
        }
    }
    
}
int main() {
    char flag = 'y';
    while(flag=='y'){
        runLength.clear();
        Timer myTimer;
        // set the bufferSize
        std::cout << "Please input the buffer's size\n";
        int bufferSize = 0;
        std::cin >> bufferSize;
        //set the losertree's size
        std::cout << "Please input the loserTree's size \n";
        int treeSize = 0;
        std::cin >> treeSize;
        std::cout << "Please input the merge order k\n";
        int k;
        std::cin >> k;
        std::cout << "Please input the mergeOrder's Size\n";
        int mBufferSize;
        std::cin >> mBufferSize;
        myTimer.startTimer();
        startSort(bufferSize,treeSize);
        mergeRuns(k, mBufferSize);
        myTimer.endTimer();
        myTimer.calculateTime();
        std::cout << "Do you want to continue?(y/n)\n";
        flag = 'n';
        std::cin >> flag;
        system("cls");
    }

    system("pause");
    return 0;
}