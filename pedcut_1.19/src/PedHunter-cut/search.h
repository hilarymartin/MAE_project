#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

/*********************************************************************/
/***     THINGS TO CHANGE FOR DIFFERENT DATABASES                  ***/
/*********************************************************************/

/* Number of people in database: */
#if !defined(POPULATION)
#define POPULATION 500000
#endif

/* Name of files with tables */
#define PERSON_TABLE_NAME   "indiv.dat"
#define CHILD_PARENT_TABLE_NAME  "relat.dat"
#define ID_TABLE_NAME  "id_table.dat"
#define GENERATION_TABLE_NAME  "gen.dat"

/* Name of pedigree for printing out LINKAGE format files */
#define PEDIGREE_ID "test_ped"

/* Column delimiter used in files where delimiter is a character */
#define DELIMITER ','
/*********************************************************************/
/***     SHOULD NOT HAVE TO CHANGE ANYTHING MORE                   ***/
/*********************************************************************/


#define MAXLINELEN 50000
#define WHITE_SPACE  " \n\r\t\v"
#define NO_RESULT   -1
#define LENGTH_CUTOFF 10  /* length of list to print on a line */

#define MALE            0
#define FEMALE          1
#define UNKNOWN_GENDER  2

#define FATHER 0
#define MOTHER 1

/* Columns for TABLES */
/* relationship table */
#define FATHER_COLUMN 1
#define MOTHER_COLUMN 2
#define CHILDREN_COLUMN 4

/* generation table */
#define PID_COLUMN 1 /*also used for person table*/
#define MAX_GENERATION_COLUMN 2

/* person table */
#define NAME_COLUMN 2
#define BIRTH_DATE_COLUMN 3
#define DEATH_DATE_COLUMN 4
#define GENDER_COLUMN 6
#define INFO_COLUMN 9

/* TRUE and FALSE may be already defined, */
/* let's not get in the way of other definitions. */
 
#if     defined(TRUE)
#undef TRUE 
#endif
 
#if     defined(FALSE)
#undef FALSE
#endif
 
typedef enum {FALSE = 0, TRUE = 1} bool; /* boolean type */
 
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif
 
#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif
 
#if !defined(NIL)
#define NIL '\0'
#endif

#define IN_GRAPH  'P'
#define OUT_GRAPH 'N'

/* Command line arguments */
int        PERSON;     /* Get some information about person */
int        PERSON2;    /* Get some information about person and person2 */
char       *INFILE_NAME;  /* input file with parameters */
char       *OUTFILE_NAME;  /* output file */
int        KEY;    /* Get key for using both parents or only males or females*/

/* store parent information for whole population */
int  **parent_information;  
int  **couple_information;  

bool PARENT_TABLE;
bool PERSON_TABLE;
char error_message[MAXLINELEN+1];

typedef struct _pid_list {
   int              pid;
   struct _pid_list *next;
} pid_list;

typedef struct _up_tree {
   int                self_pid;
   struct _up_tree    *father;
   struct _up_tree    *mother;
   struct _up_tree    *child;
} up_tree;

typedef struct _up_tree_list {
   up_tree       *node;
   struct _up_tree_list *next;
} up_tree_list;

typedef struct {
    char   *name;
    char   gender;
    char   *birth_date;
    char   *death_date;
    char   *useful_information;
} person_record;


/* store parent information for whole population */
person_record *person_information;



/* Functions in dbaccess.c for C/SQL interface and parsing arguments */
void parse_arguments(int, char **);
void suffix(char **, char *, int);
void int_suffix(int *, char *, int);
bool correct_argument(int);
void print_exit(char *);

/* Functions in families.c */
pid_list *lookup_siblings(int);
pid_list *half_siblings(int);
pid_list *uncles_and_aunts(int);
pid_list *first_cousins(int);

/* Functions in list_functions.c */
void throw_if_null(void *, char *);
int length_list(pid_list *);
void print_list(pid_list *);
void free_list(pid_list *);
pid_list *subtract_list(pid_list *, pid_list *);
void delete_pid(pid_list **, int);
void add_pid(int, pid_list **);
void append_pid(int, pid_list **);
pid_list *duplicate_list(pid_list *);
void reverse_list(pid_list **);
bool member(int, pid_list *);
void union_list(pid_list **, pid_list *);
pid_list *intersection_list(pid_list *, pid_list *);
pid_list *read_pids(char *);
void list_for_children_pids(char *, pid_list **);

/* Functions in pedigree.c */
pid_list *ancestors(int, int);
pid_list *all_common_ancestors(int, int);
pid_list *lca(int, int);
pid_list *delete_ancestors_in_list(pid_list *);
pid_list *minimal_ancestors(pid_list *);
pid_list *ancestors_length(int, int);

/* Functions for table_download */
void get_parent_information();  /* download child_parent_table */
int lookup_father(int);
int lookup_mother(int);
pid_list *lookup_spouses(int);
pid_list *lookup_spouses_with_children(int);
pid_list *lookup_children(int);
pid_list *lookup_children_for_couple(int,int);
void get_person_information();  /* download person_table */

/* Functions for printing.c */
void print_gender(int);
void print_name_id(int);
char *get_id(int);
char *get_name(int);
void print_person_information(int);
void print_rel_linkage(pid_list *, char *);
void print_list_linkage(pid_list *);
