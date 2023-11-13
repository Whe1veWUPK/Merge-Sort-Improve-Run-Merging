#include<chrono>
#include<iostream>
#include<thread>
#include<vector>
#include<queue>
#include<deque>
#include<fstream>
#include "LoserTree.hpp"
#include "FileOperator.hpp"
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
void readRunData(Buffer*freeBuffer,FileOperator*fileOperator,std::string filePath,std::deque<int>*workingBufferQueue,int runNum,int*runCurPosition,int*runStartPosition,bool isFirst,bool*isOver,int*curRunNum,int dataSize,std::queue<int>*freeIndex,int bufferSize){
    //the process has been done in intial part
    if(isFirst){
        return;
    }
    //find the run that will exhaust first
    int key = -1;
    int minValue = (((1 << 30) - 1) * 2 + 1);
    for (int i = 0; i < runNum;++i){
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
    if(key==runNum-1){
            fileOperator->writeToInputBuffer(filePath, bufferSize, runCurPosition[key], &freeBuffer[workingBufferQueue[key].back()]);
            runCurPosition[key] += freeBuffer[workingBufferQueue[key].back()].getCurSize();
    }
    else if(runCurPosition[key]+bufferSize<=runStartPosition[key+1]){
        fileOperator->writeToInputBuffer(filePath, bufferSize, runCurPosition[key], &freeBuffer[workingBufferQueue[key].back()]);
        runCurPosition[key] += bufferSize;            
    }
    else if((runCurPosition[key]+bufferSize>runStartPosition[key+1])&&(runCurPosition[key]<runStartPosition[key+1])){
        fileOperator->writeToInputBuffer(filePath, runStartPosition[key + 1] - runCurPosition[key],runCurPosition[key],&freeBuffer[workingBufferQueue[key].back()]);
        runCurPosition[key] = runStartPosition[key + 1];
    }
    else if(runCurPosition[key]==runStartPosition[key+1]){
        //the run is over
        return;
    }
    //set the over info and curRunNum info
    if(key==runNum-1){
        //if the last run is over, then return
        if(isOver[key]){
            return;
        }
        //if the last run's data has all been read, set the run's status to over 
        if(runCurPosition[key]>=dataSize&&(!isOver[key])){
            isOver[key] = true;
            --*curRunNum;
        }

    }
    else if(runCurPosition[key]==runStartPosition[key+1]){
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
    if(freeBuffer[workingBufferQueue[key].back()].getCurSize()==0){
       //the run is over
       //1.resize the buffer
        freeBuffer[workingBufferQueue[key].back()].resize(bufferSize);
       //2. pop_back it
        workingBufferQueue[key].pop_back();
       //3  push it to the freeIndex
        freeIndex->push(newKey);
    }
    else{
     //std::cout << "Update Run " << key <<" at FreeBuffer: " <<workingBufferQueue[key].back()<<"\n";
     if(key!=runNum-1){
     //std::cout << "CurPos: " << runCurPosition[key] <<" End At: " <<runStartPosition[key+1]<<"\n";
     }
    //  for (int i = 0; i < freeBuffer[newKey].getCurSize();++i){
    //     std::cout << freeBuffer[newKey].buffer[i] << " ";
    //  }
    //  std::cout << "\n";
    }
}
//k-way-merge
void startMerge(LoserTree*loserTree,Buffer*outputBuffer,Buffer*freeBuffer,std::deque<int>*workingBufferQueue,int*runNum,int*runNumNotMerged,bool*isOver,std::queue<int>*freeIndex){
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
void startSort() {
    FileOperator file;
    int dataNum = file.getDataSize("Input.txt");
    //output the data size
    std::cout << "The data size is " << dataNum << "\n";
    // set the bufferSize
    std::cout << "Please input the buffer's size\n";
    int bufferSize = 0;
    std::cin >> bufferSize;
    Buffer* buffers = new Buffer[3];
    for (int i = 0; i < 3; ++i) {
        buffers[i] = Buffer(bufferSize);
    }
    //set the losertree's size
    std::cout << "Please input the loserTree's size \n";
    int treeSize = 0;
    std::cin >> treeSize;
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
void mergeRuns(){
    //the input and output file path
    std::string inputFilePath = "Input.txt";
    std::string outputFilePath = "Output.txt";
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
    // input the k(Merge Order) and the buffer's size
    std::cout << "Please input the k(Merge Order) to continue:\n";
    int k;
    std::cin >> k;
    std::cout << "Please input the buffer's size:\n";
    int bufferSize;
    std::cin >> bufferSize;
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
    std::deque<int>*workingBufferQueue=new std::deque<int>[runNum];
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
       runCurPosition[i] = runStartPosition[i];
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
    if (runNum < k){
      //if num of run is less than k
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
          if ((i != runNum - 1) && (runCurPosition[i] + bufferSize <= runStartPosition[i + 1])){
              f.writeToInputBuffer(outputFilePath, bufferSize, runCurPosition[i], &freeBuffer[workingBufferQueue[i].front()]);
              runCurPosition[i] += bufferSize;
            }
          else if (i == runNum - 1){
              f.writeToInputBuffer(outputFilePath, bufferSize, runCurPosition[i], &freeBuffer[workingBufferQueue[i].front()]);
              runCurPosition[i] += bufferSize;
            }
          else if (runCurPosition[i] + bufferSize > runStartPosition[i + 1]){
              f.writeToInputBuffer(outputFilePath, runStartPosition[i + 1] - runCurPosition[i], runCurPosition[i], &freeBuffer[workingBufferQueue[i].front()]);
              runCurPosition[i] = runStartPosition[i + 1];
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
       if(key==runNum-1){
            f.writeToInputBuffer(outputFilePath, bufferSize, runCurPosition[key], &freeBuffer[index]);
            runCurPosition[key] += bufferSize;
       }
       else if(runCurPosition[key]+bufferSize<=runStartPosition[key+1]){
            f.writeToInputBuffer(outputFilePath, bufferSize, runCurPosition[key], &freeBuffer[index]);
            runCurPosition[key] += bufferSize;            
       }
       else if(runCurPosition[key]+bufferSize>runStartPosition[key+1]){
            f.writeToInputBuffer(outputFilePath, runStartPosition[key + 1] - runCurPosition[key],runCurPosition[key],&freeBuffer[index]);
            runCurPosition[key] = runStartPosition[key + 1];
       }
       //set activeOutputBuffer
       int activeOutputBuffer = 0;
       int curNum = runNum;
       bool isFirst = true;
       bool test = true;
    //    for (int i = 0; i < runNum;++i){
    //         std::cout << "workingBuffer " << (i) << "'s size is: " << workingBufferQueue[i].size() << "\n";
    //    }
       //startMerge
       std::fstream fout;
       if (isFirst) {
        //如果是某一轮次的第一个 outputBuffer输出 则重写文档
        fout.open("Input.txt", std::ios::out | std::ios::trunc);
       }
       else {
        
        fout.open("Input.txt", std::ios::out | std::ios::app);
       }
       fout.close();
       int totalRunNum = runNum;
       while((runNotFinishedMergeNum>0)){
        // for (int i = 0; i < totalRunNum;++i){
        //       std::cout << "Run " << i << "'s data: \n";
        //       for (int j = 0; j <workingBufferQueue[i].size();++j){
        //             std::cout << "FreeBuffer " << workingBufferQueue[i].at(j) <<"\n";
        //             for (int k = 0; k < freeBuffer[workingBufferQueue[i].at(j)].getCurSize();++k){
        //                 std::cout << freeBuffer[workingBufferQueue[i].at(j)].buffer[k] << " ";
        //             }
        //             std::cout << "\n";
        //       }
        //       std::cout << "\n";
        // }
        
        // thread 1: do the k-way-merge in memory
        std::thread t1(startMerge, &loserTree, &outputBuffer[activeOutputBuffer], freeBuffer, workingBufferQueue, &runNum,&runNotFinishedMergeNum,isOver, &freeIndex);
        // thread 2: write the another output buffer's data to the disk
        std::thread t2(writeToDisk, &outputBuffer[(activeOutputBuffer + 1) % 2], &f, isFirst, inputFilePath);
        // wait thread 1 to finish
        t1.join();
        // thread 3: read the run's data that it will first exhaust
        std::thread t3(readRunData, freeBuffer, &f, outputFilePath, workingBufferQueue, totalRunNum, runCurPosition, runStartPosition, isFirst,isOver,&runNum,dataSize,&freeIndex, bufferSize);
        // wait thread 2 and thread 3 to finish
        t2.join();
        t3.join();
        isFirst = false;
        // std::cout << "Output buffer generated!\n";
        // test = false;
        activeOutputBuffer = (activeOutputBuffer + 1) % 2;
        }
       // write the last output buffer into the disk
        writeToDisk(&outputBuffer[(activeOutputBuffer + 1) % 2], &f, false, inputFilePath);

        

  }
    
}
int main() {
    char flag = 'y';
    while(flag=='y'){
        runLength.clear();
        startSort();
        //getBestMergeSequence();
        mergeRuns();
        std::cout << "Do you want to continue y/n(Yes/No)?\n";
        std::cin >> flag;
        system("cls");
    }

    system("pause");
    return 0;
}