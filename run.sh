echo -------UNIT TEST --- OPENMP---------
make openmp; ./forth.out 7000 1 ; ./forth.out 7000 4 ; ./forth.out 7000 6;
make clean;
echo -------UNIT TEST --- PTHREADS ------
make pthread; ./pthread.out 7000 1 ; ./pthread.out 7000 4; ./pthread.out 7000 6;
make clean;
echo -------------EXITING----------------
