#include <iostream>
#include <vector>
#include <cstdlib>
#include <time.h>
#include <mpi.h>

/*
Some characterstics of the program: 
1. Each process determines their part of the domain
   [Min] 0----4 5----9 10----14 15----20 [Max]
           (p1)   (p2)    (p3)    (p4)
2. Each process initializes exactly N walker, all 
   start at the first value of their domain
3. Each walker has two associated integer values:
   - current position of the walker
   - numner of steps left to take
4. Walker starts traversing through the domain and 
   passed to other processes until they have completed
   their walk
5. The process terminitaes when all walkers have finished
*/

using namespace std;

typedef struct {
  int location;
  int num_steps_left_in_walk;
} Walker;

/*
  Domain Decomposition: it takes the domian size and find the
  appropriate subdomain for MPI process
*/
void decompose_domain(int domain_size, int world_rank,
                      int world_size, int* subdomain_start,
                      int* subdomain_size) {
  if (world_size > domain_size) {
    // Don't worry about this special case. Assume the domain size
    // is greater than the world size.
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  *subdomain_start = domain_size / world_size * world_rank;
  *subdomain_size = domain_size / world_size;
  if (world_rank == world_size - 1) {
    // Give remainder to last process
    *subdomain_size += domain_size % world_size;
  }
}

/*
 Our initialization function, called initialize_walkers, takes the
 subdomain bounds and adds walkers to an incoming_walkers vector. 
*/
void initialize_walkers(int num_walkers_per_proc, int max_walk_size,
                        int subdomain_start,
                        vector<Walker>* incoming_walkers) {
  Walker walker;
  for (int i = 0; i < num_walkers_per_proc; i++) {
    // Initialize walkers at the start of the subdomain
    walker.location = subdomain_start;
    walker.num_steps_left_in_walk =
      (rand() / (float)RAND_MAX) * max_walk_size;
    incoming_walkers->push_back(walker);
  }
}

/*
  This function is responsible for progressing the walker until
  it has finished its walk. If it goes out of local bounds, it
  is added to the outgoing_walkers vector.
*/
void walk(Walker* walker, int subdomain_start, int subdomain_size,
          int domain_size, vector<Walker>* outgoing_walkers) {
  while (walker->num_steps_left_in_walk > 0) {
    if (walker->location == subdomain_start + subdomain_size) {
      // Take care of the case when the walker is at the end
      // of the domain by wrapping it around to the beginning
      if (walker->location == domain_size) {
        walker->location = 0;
      }
      outgoing_walkers->push_back(*walker);
      break;
    } else {
      walker->num_steps_left_in_walk--;
      walker->location++;
    }
  }
}

void send_outgoing_walkers(vector<Walker>* outgoing_walkers,
                           int world_rank, int world_size) {
  // Send the data as an array of MPI_BYTEs to the next process.
  // The last process sends to process zero.
  MPI_Send((void*)outgoing_walkers->data(),
           outgoing_walkers->size() * sizeof(Walker), MPI_BYTE,
           (world_rank + 1) % world_size, 0, MPI_COMM_WORLD);
  // Clear the outgoing walkers list
  outgoing_walkers->clear();
}

void receive_incoming_walkers(vector<Walker>* incoming_walkers,
                              int world_rank, int world_size) {
  // Probe for new incoming walkers
  MPI_Status status;
  // Receive from the process before you. If you are process zero,
  // receive from the last process
  int incoming_rank =
    (world_rank == 0) ? world_size - 1 : world_rank - 1;
  MPI_Probe(incoming_rank, 0, MPI_COMM_WORLD, &status);
  // Resize your incoming walker buffer based on how much data is
  // being received
  int incoming_walkers_size;
  MPI_Get_count(&status, MPI_BYTE, &incoming_walkers_size);
  incoming_walkers->resize(incoming_walkers_size / sizeof(Walker));
  MPI_Recv((void*)incoming_walkers->data(), incoming_walkers_size,
           MPI_BYTE, incoming_rank, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
}


/*
We have to tie all these function together as follows:
    1. Initialize the walkers.
    2. Progress the walkers with the walk function.
    3. Send out any walkers in the outgoing_walkers vector.
    4. Receive new walkers and put them in the incoming_walkers vector.
    5. Repeat steps two through four until all walkers have finished.
*/
int main(int argc, char** argv) {
  int domain_size;
  int max_walk_size;
  int num_walkers_per_proc;

  if (argc < 4) {
    cerr << "Usage: random_walk domain_size max_walk_size "
         << "num_walkers_per_proc" << endl;
    exit(1);
  }
  domain_size = atoi(argv[1]);
  max_walk_size = atoi(argv[2]);
  num_walkers_per_proc = atoi(argv[3]);

  MPI_Init(NULL, NULL);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  srand(time(NULL) * world_rank);
  int subdomain_start, subdomain_size;
  vector<Walker> incoming_walkers, outgoing_walkers;

  // Find your part of the domain
  decompose_domain(domain_size, world_rank, world_size,
                   &subdomain_start, &subdomain_size);
  // Initialize walkers in your subdomain
  initialize_walkers(num_walkers_per_proc, max_walk_size, subdomain_start,
                     &incoming_walkers);

  cout << "Process " << world_rank << " initiated " << num_walkers_per_proc
       << " walkers in subdomain " << subdomain_start << " - "
       << subdomain_start + subdomain_size - 1 << endl;

  // Determine the maximum amount of sends and receives needed to
  // complete all walkers
  int maximum_sends_recvs = max_walk_size / (domain_size / world_size) + 1;
  for (int m = 0; m < maximum_sends_recvs; m++) {
    // Process all incoming walkers
    for (int i = 0; i < incoming_walkers.size(); i++) {
       walk(&incoming_walkers[i], subdomain_start, subdomain_size,
            domain_size, &outgoing_walkers);
    }
    cout << "Process " << world_rank << " sending " << outgoing_walkers.size()
         << " outgoing walkers to process " << (world_rank + 1) % world_size
         << endl;
    if (world_rank % 2 == 0) {
      // Send all outgoing walkers to the next process.
      send_outgoing_walkers(&outgoing_walkers, world_rank,
                            world_size);
      // Receive all the new incoming walkers
      receive_incoming_walkers(&incoming_walkers, world_rank,
                               world_size);
    } else {
      // Receive all the new incoming walkers
      receive_incoming_walkers(&incoming_walkers, world_rank,
                               world_size);
      // Send all outgoing walkers to the next process.
      send_outgoing_walkers(&outgoing_walkers, world_rank,
                            world_size);
    }
    cout << "Process " << world_rank << " received " << incoming_walkers.size()
         << " incoming walkers" << endl;
  }
  cout << "Process " << world_rank << " done" << endl;
  MPI_Finalize();
  return 0;
}
