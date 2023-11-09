#include<chrono>
#include<iostream>
#include<thread>
#include<vector>
#include<queue>
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
//Thread 1: read data from the disk
void readFromDisk(Buffer* inputBuffer, FileOperator* fileOperator, int startPosition, std::string filePath) {
    //get inputBuffer's Size 
    int predictSize = inputBuffer->getBufferSize();
    //
    inputBuffer->resize(predictSize);
    //read info from file
    fileOperator->writeToInputBuffer(filePath, predictSize, startPosition, inputBuffer);
}
//Thread 2: generate Run
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
//Thread 3: write data to the disk
void writeToDisk(Buffer* outputBuffer, FileOperator* FileOperator, bool isFirst, std::string  filePath) {

    //get outputBuffer's size
    int predictSize = outputBuffer->getCurSize();
    if (outputBuffer->Empty()) {
        std::cerr << "The outputBuffer is empty! \n";
        return;
    }
    //write the buffer's data to disk
    FileOperator->writeToFile(filePath, predictSize, isFirst, outputBuffer);
    if (isFirst) {
        isFirst = false;
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
    std::cout << "initial buffer already\n";
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

    FileOperator f;
    //get the total data size
    int dataSize = 0;
    for (int i = 0; i < runLength.size();++i){
        dataSize += runLength[i];
    }
    int runNum = runLength.size();
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
    //flag to show if the buffer is free
    bool *isFree = new bool[2 * k];
    //cur pin
    int curPin = 0;
    for (int i = 0; i < 2 * k;++i){
        freeBuffer[i] = Buffer(bufferSize);
        isFree[i] = true;
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
    
    Buffer**workingBuffer=new Buffer*[k];
    for (int i = 0; i < k;++i){
        workingBuffer[i] = new Buffer[2];
    }
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
      //initiate working buffer
       for (int i = 0; i < runNum;++i){
         workingBuffer[i][0] = freeBuffer[i];
         isFree[i] = false;
       }
      //read data from disk to working buffer
       for (int i = 0; i < runNum;++i){
          if((i!=runNum-1)&&(runCurPosition[i]+bufferSize<=runStartPosition[i+1])){
                f.writeToInputBuffer(outputFilePath, bufferSize, runCurPosition[i], &workingBuffer[i][0]);
                runCurPosition[i] += bufferSize;
          }
          else if(i==runNum-1){
                f.writeToInputBuffer(outputFilePath, bufferSize, runCurPosition[i], &workingBuffer[i][0]);
                runCurPosition[i] += bufferSize;
          }
          else if(runCurPosition[i]+bufferSize>runStartPosition[i+1]){
                f.writeToInputBuffer(outputFilePath, runStartPosition[i + 1] - runCurPosition[i], runCurPosition[i], &workingBuffer[i][0]);
                runCurPosition[i] = runStartPosition[i + 1];
          }
       }
       //initiate data
       int *data = new int[runNum];
       for (int i = 0; i < runNum;++i){
          data[i] = workingBuffer[i][0].buffer[0];
       }

       for (int i = 0; i < runNum;++i){
          std::cout << data[i] << " ";
       }
       std::cout << "\n";
       LoserTree loserTree(data, runNum);
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