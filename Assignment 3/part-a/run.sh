g++ -std=c++11 -I include/detail -I include/ -o mr-pr-cpp.o  mr-pr-cpp.cpp -lboost_system -lpthread -lboost_iostreams -lboost_filesystem;
./mr-pr-cpp.o $1.txt -o $1-pr-cpp.txt
