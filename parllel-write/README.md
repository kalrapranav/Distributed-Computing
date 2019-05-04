#Parllel File Writing

  The following program allows multiple processors to write on a same binary file using MPI-IO

##Working:
    The follwing program write data to a binary file via MPI I/O. I need process 0 to write a short header, then I need the whole range of processes to write their own pieces of the array indicated by the header. Then I need process 0 to write another header, followed by all processes writing their pieces of the next array, etc. 

 
    Source: 
    [href="https://stackoverflow.com/questions/37838228/mpi-i-o-mix-of-single-and-multiple-process-output] (href="https://stackoverflow.com/questions/37838228/mpi-i-o-mix-of-single-and-multiple-process-output)
  
    Some suggestions for more efficiency:
      * You can consult the status object for how many bytes were written, instead of getting the position and translating into bytes.
      * If you have the memory to hold all the data before you write, you could describe your I/O with an MPI datatype (admittedly, one that might end up being a pain to create). Then all processes would issue a single collective call.
      * You should use collective I/O instead of independent I/O. A "quality library" should be able to give you equal if not better performance (and if not, you could raise the issue with your MPI implementation).
      * If the processes have different amounts of data to write, MPI_EXSCAN is a good way to collect who has what data. Then you can call MPI_FILE_WRITE_AT_ALL to the correct offset in the file.

![parllel-write](https://user-images.githubusercontent.com/19777060/57171544-1a164800-6dca-11e9-94bf-4e0c1b2f96b0.jpg)