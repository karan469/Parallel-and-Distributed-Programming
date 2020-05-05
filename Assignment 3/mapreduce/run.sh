g++ -g -std=c++11 -I include/detail -I include/ -o "$1".o  "$1".cpp table.cpp -lboost_system -lpthread -lboost_iostreams -lboost_filesystem
