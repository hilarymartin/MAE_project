#include "search.h"

typedef unsigned int HashKey;

typedef struct {
  HashKey key;
  int p1;
  int p2;
  double kinship_value;
} KinPair;

typedef struct {
  KinPair *slots;
  int tablesize;
  int tablebits;
  int slots_used;
} KinshipTable;

#define NOFREESLOTS ((unsigned int) -1)
#define NO_INFO ((unsigned int) -1)
#define NO_VALUE -1.0

#define TABLESIZE  (POPULATION * 4)

int *generation_information;

KinshipTable *NewKinshipTable(int);
bool KinshipTableInsert(KinshipTable *, int, int, double);
double KinshipTableLookup(KinshipTable *, int, int);

double inbreeding(int, KinshipTable *);
double kinship(int,int, KinshipTable *);

void get_generation_information();
char *get_token(char *,int);
