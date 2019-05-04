[Output with 2 processes: mpirun -n 2 ./ring]
Process 1 received token -1 from process 0
Process 0 received token -1 from process 1

[Output with 4 processes: mpirun -n 4 ./ring]
Process 1 received token -1 from process 0
Process 2 received token -1 from process 1
Process 3 received token -1 from process 2
Process 0 received token -1 from process 3

[Output with 10 processes: mpirun -n 10 ./ring]
Process 1 received token -1 from process 0
Process 2 received token -1 from process 1
Process 3 received token -1 from process 2
Process 4 received token -1 from process 3
Process 5 received token -1 from process 4
Process 6 received token -1 from process 5
Process 7 received token -1 from process 6
Process 0 received token -1 from process 9
Process 8 received token -1 from process 7
Process 9 received token -1 from process 8
