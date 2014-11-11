/*
 * asp.c
 *
 * Usage: asp -F<infile> -K<parents>
 *
 * Find all shortest paths pedigree for persons in <infile>
 * -K determines Key with 1 for males only, 2 for females only
 * and 0 for both male and female ancestors for asp.
*/

#include "search.h"

#define REL_FIELDS 3    /*Number of fields of record needed for relationships*/
#define SELF 0
#define FATHER_INDEX 1
#define MOTHER_INDEX 2

/* different printout formats for ASP pedigree */
#define SHOW_PATHS 1
#define FOR_MAKEPED 2

/*structure for node in a tree*/
typedef struct _tree_node {
    int  self;
    struct _tree_node *father;
    struct _tree_node *mother;
    struct _tree_node_list *children;
} tree_node;

/*structure for list of tree_nodes*/
typedef struct _tree_node_list {
    tree_node *this;
    struct _tree_node_list *next;
} tree_node_list;

static tree_node *short_ancestors_tree(int, pid_list **);
static tree_node *find_tree_node_in_list(int, tree_node_list *);
static tree_node *find_tree_node_in_uptree(int, tree_node *);
static pid_list *list_from_tree(tree_node *);
static void merge_tree(tree_node *, tree_node *);
static tree_node *get_tree_node(int);
static tree_node_list *get_tree_node_list(tree_node *);
static void add_child(tree_node *, tree_node *);
static void add_node_to_tree_list(tree_node *, tree_node_list **);
static void free_tree_node_list(tree_node_list *);
static void print_tree_linkage(tree_node *, int);
static void add_child_parent_link(int **, tree_node *, int);
static void add_link(int, int, int **, int, int);
static int lookup_gender(int, int);
static pid_list *list_individuals_down(tree_node *);

main(argc, argv)
int             argc;
char            *argv[];
{
    pid_list *given_people = NULL;  /*list of pids for people specified to be 
                                        in asp pedigree*/
    int number_of_given_people;      /*number of people specified*/
    pid_list *temp = NULL;       /*pointer into list of given people*/
    int index;                   /*counter for list of given people*/
    tree_node **terminals = NULL;  /*an array of ancestor trees for
                                     given people (terminals), one per person*/
    pid_list *common_ancestors = NULL; /*list of common ancestors*/
    pid_list *result = NULL;  /*list of common ancestors
                                with non-"least" ancestors pruned*/
    pid_list *ancestors_list = NULL;  /*ancestors for one person*/
    pid_list *plist = NULL;  /*placeholder list for partial list of
                               common ancestors*/
    tree_node *leaf = NULL;  /*tree node storing a terminal*/
    int root;        /* id of root node in a tree of some ancestors*/
    tree_node *root_node = NULL;  /*tree_node for the root of a tree*/
    tree_node *asp_graph = NULL;  /* graph of all shortest paths*/

    parse_arguments(argc, argv);

    if (!(correct_argument('F') && correct_argument('K')))
        print_exit("Usage: asp -F<infile> -K<Key>\n\tKey: 0 for both, 1 for males only, 2 for females only\n");

    get_person_information();
    get_parent_information();

    given_people = read_pids(INFILE_NAME);
    if (given_people == NULL)
    {
        fprintf(stderr,"Specify people for finding ASP in the input file\n");
        exit(EXIT_FAILURE);
    }
 
    /* Count number of individuals in the input file */
    number_of_given_people = 0;
    for (temp = given_people; temp != NULL; temp = temp->next)
        number_of_given_people++;

    /* allocate space for nodes for given people */
    terminals = (tree_node **) malloc(number_of_given_people *
				sizeof(tree_node *));
    if (terminals == NULL)
    {
        fprintf(stderr,"Cannot allocate space for nodes for given people\n");
        exit(EXIT_FAILURE);
    }
 
    /* find tree of ancestors and common ancestors */
    common_ancestors = NULL;

    printf("Searching for ancestors for %d ...\n", given_people->pid);

    /* Find ancestor tree for first given person and
       Initialize common ancestors with all ancestors for first given person */

    terminals[0] = short_ancestors_tree(given_people->pid, &common_ancestors);

    /* Find ancestor tree for remaining given person and
       find common ancestors as intersection of all ancestors for each.
       length of given_people list = number_of_given_people so
       temp and index go together */

    for (temp = given_people->next, index = 1;
         ((temp != NULL) && (common_ancestors != NULL));
         temp = temp->next, index++)
    {
        /* ancestors for given person being considered */
        ancestors_list = NULL;
        printf("Searching for ancestors for %d ...\n", temp->pid);
        terminals[index] = short_ancestors_tree(temp->pid, &ancestors_list);

        /* introducing plist so that old common_ancestors list can be freed */
        plist = common_ancestors;
        common_ancestors = intersection_list(plist, ancestors_list);
        free_list(plist);
        free_list(ancestors_list);
    }

    if (common_ancestors == NULL)
    {
        printf("Can't find a common ancestor for ");
        print_list(given_people);
    }
    else
    {
        printf("All common ancestors for ");
        print_list(given_people);
        printf("are ");
        print_list(common_ancestors);

        /* Find minimal from list of all ancestors */
        result = duplicate_list(common_ancestors);
        for (temp = common_ancestors; temp != NULL; temp = temp->next)
        {
            /* Since common ones are common to all, we can use any ancestor 
               tree, in particular, terminals[0].
               terminal[0] is a leaf, find the common ancestor in tree
               above this leaf */
            leaf = find_tree_node_in_uptree(temp->pid, terminals[0]);

            /* find list of ancestors for leaf */
            ancestors_list = list_from_tree(leaf);
            delete_pid(&ancestors_list, leaf->self);
        
            /* remove any in common with all ancestors */
            for (plist = ancestors_list; plist != NULL; plist = plist->next)
                delete_pid(&result, plist->pid);
 
            free_list(ancestors_list);
        }
        free_list(common_ancestors);
        common_ancestors = result;

        if (common_ancestors == NULL)
        {
            sprintf(error_message,"Cannot have common but no minimal\n");
            print_exit(error_message);
        }

        if (common_ancestors->next == NULL)
            printf("Minimal ancestor that will serve as root is %d",
		    common_ancestors->pid);
        else
        {
            printf("Minimal ancestors that will serve as roots are ");
            print_list(common_ancestors);
        }
 
        /* Merge trees above each given person for each common ancestor */
        for (temp = common_ancestors; temp != NULL; temp = temp->next)
        {
            root = temp->pid;
            asp_graph = get_tree_node(root);

            for (index = 0; index < number_of_given_people; index++)
            {
                root_node = find_tree_node_in_uptree(root, terminals[index]);
                merge_tree(root_node, asp_graph);
            }
            printf("\n");
            printf("\t Pedigree showing shortest paths to common ancestor %d\n", root);
            printf("\n");
            print_tree_linkage(asp_graph, SHOW_PATHS);

            printf("\n");
            printf("\t Pedigree to common ancestor %d for linkage studies\n", root);
            printf("\n");
            print_tree_linkage(asp_graph, FOR_MAKEPED);
        }
    }
}

/* Return graph of ancestors for leaf and 
   list of all ancestors in ancestors_list */
tree_node *short_ancestors_tree(leaf, ancestors_list)
int leaf;
pid_list **ancestors_list;
{
    tree_node *tree = NULL; /*tree node representing the initial leaf*/
    pid_list *individuals_in_tree = NULL; /*list of pids of individuals in tree*/
    pid_list *next_list = NULL;  /*next list of parents of expanded nodes*/
    tree_node_list *nodes_to_expand = NULL; /*set of nodes whose parents
                                              may not have been added yet*/
    tree_node_list *next_set = NULL;        /*next set of nodes to expand*/
    tree_node_list *temp;      /*pointer to next node to expand in list*/
    tree_node *child_node;      /*node to expand as a tree_node*/
    int pid;                   /*placeholder for a new found pid*/

    add_pid(leaf, &individuals_in_tree);
    tree = get_tree_node(leaf);
    nodes_to_expand = get_tree_node_list(tree);

    while (nodes_to_expand != NULL)
    {
        next_set = NULL;  /* next set of nodes to expand */
        next_list = NULL;  /* list of parents of nodes_to_expand */

        /* find parents for each node */
        for (temp = nodes_to_expand; temp != NULL; temp = temp->next)
        {
            child_node = temp->this;
    
            if ((KEY == 0) || (KEY == 1))
                pid = lookup_father(child_node->self);
            else
                pid = NO_RESULT;
             
            if (pid <= NO_RESULT)
                child_node->father = NULL;
            else
            {
                /* if already in list, the path we found is longer.
                   We won't add link from parent to child as we want
                   shortest paths from parent -> child but need to be
		   able to find all ancestors */
                if (member(pid, individuals_in_tree))
                   child_node->father = find_tree_node_in_uptree(pid, tree);
                else
                {
                   if (member(pid, next_list)) /* another same length path */
                   {
                      child_node->father = find_tree_node_in_list(pid,next_set);
                      add_child(child_node->father, child_node);
                   }
                   else /* introduce node for father */
                   {
                       child_node->father = get_tree_node(pid);
                       add_child(child_node->father, child_node);
                       add_node_to_tree_list(child_node->father, &next_set);
                       add_pid(pid, &next_list);
                   } 
                }
            }

            /* Do same thing for mother */
            if ((KEY == 0) || (KEY == 2))
                pid = lookup_mother(child_node->self);
            else
                pid = NO_RESULT;
             
            if (pid <= NO_RESULT)
                child_node->mother = NULL;
            else
            {
                if (member(pid, individuals_in_tree))
                   child_node->mother = find_tree_node_in_uptree(pid, tree);
                else
                {
                   if (member(pid, next_list))
                   {
                      child_node->mother = find_tree_node_in_list(pid,next_set);
                      add_child(child_node->mother, child_node);
                   }
                   else 
                   {
                       child_node->mother = get_tree_node(pid);
                       add_child(child_node->mother, child_node);
                       add_node_to_tree_list(child_node->mother, &next_set);
                       add_pid(pid, &next_list);
                   }
                }
            }
        }

        union_list(&individuals_in_tree,next_list);
        free_list(next_list);
        free_tree_node_list(nodes_to_expand);

        nodes_to_expand = next_set;
    }

    *ancestors_list = individuals_in_tree;
    return(tree);
}

/*Find and return the tree_node the tree_node whose
  self field matches pid from a list of tree_nodes*/
tree_node *find_tree_node_in_list(pid, list)
int pid;
tree_node_list *list;
{
    tree_node_list *temp;

    for (temp = list; temp != NULL; temp = temp->next)
        if (temp->this->self == pid)
            return(temp->this);

    return(NULL);
}

/*Find and return a tree_node whose self field matches pid and
  is an ancestor of child*/
tree_node *find_tree_node_in_uptree(pid, child)
int pid;
tree_node *child;
{
    tree_node *possible_result; /*possible result in father's ancestors*/

    if (child == NULL)
        return(NULL);

    if (child->self == pid)
        return(child);

    /* There is at most one node in tree above child for pid.
       If above father, we don't need to look at tree above mother */
    possible_result = find_tree_node_in_uptree(pid, child->father);
    if (possible_result != NULL)
        return(possible_result);
    else
        return(find_tree_node_in_uptree(pid, child->mother));
}

/* Compute and return a list of all ancestors of leaf */
pid_list *list_from_tree(leaf)
tree_node *leaf;
{
    pid_list *result = NULL; /*result list to return*/
    pid_list *one_parent_list; /*partial list used in intermediate step*/

    if (leaf == NULL)
        return(NULL);

    add_pid(leaf->self, &result);
    one_parent_list = list_from_tree(leaf->father);
    if (one_parent_list != NULL)
        union_list(&result, one_parent_list);
    one_parent_list = list_from_tree(leaf->mother);
    if (one_parent_list != NULL)
        union_list(&result, one_parent_list);

    return(result);
}

/* Merge all descendant paths originating at <from> into the tree node <to> 
comprising of result so far. from->self == to->self */
void merge_tree(from, to)
tree_node *from;
tree_node *to;
{
    tree_node_list *temp;  /*pointer into list of children of from*/
    tree_node_list *child_list; /*list of children of to*/
    tree_node_list *new; /*record for new child in a list*/
    tree_node *new_child; /*record for new child as one entry*/

    if (from == NULL)
    {
        fprintf(stderr,"Copying from NULL node?\n");
        exit(EXIT_FAILURE);
    }

    for (temp = from->children; temp != NULL; temp = temp->next)
    {
        child_list = to->children; 

        /* find the child if it is already present in to->children */
        while ((child_list != NULL) && 
               (child_list->this->self != temp->this->self))
            child_list = child_list->next;

        if (child_list == NULL)
        {
            /* Create node for child */
            new_child = get_tree_node(temp->this->self);

            /* is from a father or mother of child */
            if ((temp->this->father != NULL) &&
                (temp->this->father->self == from->self))
                new_child->father = to;
            else
            {
                if ((temp->this->mother != NULL) &&
                    (temp->this->mother->self == from->self))
                    new_child->mother = to;
                else
                {
                    fprintf(stderr,"How did you find %d\n",new_child->self);
                    exit(EXIT_FAILURE);
                }
            }

            /* Add child to list of children */
            new = get_tree_node_list(new_child);
            new->next = to->children;
            to->children = new;

            merge_tree(temp->this, to->children->this);
        }
        else
            merge_tree(temp->this, child_list->this);
    }
}

/*Allocate a new item of type tree_node and fill its
  self field with person, initialize other fields to NULL,
  return new item*/
tree_node *get_tree_node(person)
int person;
{
    tree_node *new_tree_node = NULL;
 
    new_tree_node = (tree_node *) malloc(sizeof(tree_node));
 
    if (new_tree_node == NULL)
    {
        fprintf(stderr,"Cannot allocate space for tree_node in solution\n");
        exit(EXIT_FAILURE);
    }
 
    new_tree_node->self = person;
    new_tree_node->father = NULL;
    new_tree_node->mother = NULL;
    new_tree_node->children = NULL;
 
    return(new_tree_node);
}

/*Allocate a new item of type tree_node_list and fill its
  data field with one_node, return new item*/
tree_node_list *get_tree_node_list(one_node)
tree_node *one_node;
{
    tree_node_list *new_tree_node_list = NULL; /*new item to allocate*/
 
    new_tree_node_list = (tree_node_list *) malloc(sizeof(tree_node_list));
 
    if (new_tree_node_list == NULL)
    {
        fprintf(stderr,"Cannot allocate space for tree_node_list \n");
        exit(EXIT_FAILURE);
    }
 
    new_tree_node_list->this = one_node;
    new_tree_node_list->next = NULL;
 
    return(new_tree_node_list);
}

/*Insert a new child at the beginning of parent's list of children*/
void add_child(parent, child)
tree_node *parent;
tree_node *child;
{
    tree_node_list *new_list; /*new beginning of list of children*/

    new_list = get_tree_node_list(child);
    new_list->next = parent->children;

    parent->children = new_list;
}

/*Insert a new node at the beginning of a list of tree nodes*/
void add_node_to_tree_list(node, list)
tree_node *node;
tree_node_list **list;
{
    tree_node_list *new; /*new entry to insert*/

    new = get_tree_node_list(node);
    new->next = *list;
  
    *list = new;
}

/*Deallocate all elements in a list of tree_nodes*/
void free_tree_node_list(list)
tree_node_list *list;
{
    if (list != NULL)
    {
        free_tree_node_list(list->next);
        free(list);
    }
}

/*Print out lines containing all indivduals in graph
  rooted at root in a format like LINKAGE format. If 
  the argument PRINTOUT_FORMAT == FOR_MAKEPED then the 
  output is LINKAGE format suitable for input to the makeped
  pre-processor program. If, instead, PRINTOUT_FORMAT == SHOW_PATHS,
  then each line is in LINKAGE format, but an individual may be duplicated
  once for being on paths and once as a required spouse. Identifiers
  of individuals on required paths start with P; individuals included only
  to fulfill the LINKAGE requirement that a person have either 0 parents
  or 2 parents have identifiers that start with N.*/
void print_tree_linkage(root, PRINTOUT_FORMAT)
tree_node *root;
int PRINTOUT_FORMAT;
{
    int **info;  /*array of information about people*/
    pid_list *individuals = NULL; /*list of individuals to print about*/
    pid_list *temp;  /*loop index over individuals list*/
    int index, i;  /*loop indices over people and fields*/
    int  count;  /*count of number of people*/
    int person; /*placeholder for person to print*/
    pid_list *others_father; /*list of extra parents needed to give everyone 2 parents*/
    pid_list *others_mother; /*list of extra parents needed to give everyone 2 parents*/
    int father_id, mother_id; /*placeholders for father and mother indices for
                                adding extra parents to others*/

    if (root != NULL)
    {
        if ((PRINTOUT_FORMAT != SHOW_PATHS) &&
            (PRINTOUT_FORMAT != FOR_MAKEPED))
        {
            fprintf(stderr, "Prinout format %d not supported in print_tree_linkage\n", PRINTOUT_FORMAT);
            exit(EXIT_FAILURE);
        }

        /* Find sorted list of individuals other than root */
        individuals = list_individuals_down(root);
        /*root is deleted from list because in printing its parents
          are known to be 0 0 */
        delete_pid(&individuals, root->self);

        count = 0;
        for (temp = individuals; temp != NULL; temp = temp->next)
             count++;

        /* space for storing info for self, father, mother */
        info = (int **)malloc(count * sizeof(int *));
        if (info == NULL)
        {
            fprintf(stderr,"Cannot allocate space for collecting print info\n");
            exit(EXIT_FAILURE);
        }

        for (index = 0; index < count; index++)
        {
            info[index] = (int *)malloc(REL_FIELDS * sizeof(int));
            if (info[index] == NULL)
            {
                fprintf(stderr,"Cannot allocate space for info[%d]\n",index);
                exit(EXIT_FAILURE);
            }
        }
        
        for (index = 0; index < count; index++)
            for (i = 0; i < REL_FIELDS; i++)
                info[index][i] = 0;

        /* initialize self in each */
        for (temp = individuals, index = 0; 
             ((temp != NULL) && (index < count)); 
             temp = temp->next, index++)
            info[index][SELF] = temp->pid;

        /* put parent information */
        add_child_parent_link(info, root, count);

        /* parents which are not in result graph but have to be in pedigree
           output because of LINKAGE format */
        others_father = NULL;
        others_mother = NULL;

        /* print root */
        if (PRINTOUT_FORMAT == SHOW_PATHS)
            printf("%s        P%d\t0\t0", PEDIGREE_ID, root->self);
        if (PRINTOUT_FORMAT == FOR_MAKEPED)
            printf("%s        %d\t0\t0", PEDIGREE_ID, root->self);
        print_gender(root->self);
        print_name_id(root->self);

        /* print individuals in graph */
        for (index = 0; index < count; index++)
        {
            person = info[index][0];
            if (PRINTOUT_FORMAT == SHOW_PATHS)
                printf("%s        P%d", PEDIGREE_ID, person);
            if (PRINTOUT_FORMAT == FOR_MAKEPED)
                printf("%s        %d", PEDIGREE_ID, person);

            if (info[index][FATHER_INDEX] != 0)
            {
                if (PRINTOUT_FORMAT == SHOW_PATHS)
                    printf("\tP%d", info[index][FATHER_INDEX]);
                if (PRINTOUT_FORMAT == FOR_MAKEPED)
                    printf("\t%d", info[index][FATHER_INDEX]);
            }
            else
            {
                father_id = lookup_father(person); 
                add_pid(father_id, &others_father);
                if (PRINTOUT_FORMAT == SHOW_PATHS)
                    printf("\tN%d", father_id);
                if (PRINTOUT_FORMAT == FOR_MAKEPED)
                    printf("\t%d", father_id);
            }

            if (info[index][MOTHER_INDEX] != 0)
            {
                if (PRINTOUT_FORMAT == SHOW_PATHS)
                    printf("\tP%d", info[index][MOTHER_INDEX]);
                if (PRINTOUT_FORMAT == FOR_MAKEPED)
                    printf("\t%d", info[index][MOTHER_INDEX]);
            }
            else
            {
                mother_id = lookup_mother(person); 
                add_pid(mother_id, &others_mother);
                if (PRINTOUT_FORMAT == SHOW_PATHS)
                    printf("\tN%d", mother_id);
                if (PRINTOUT_FORMAT == FOR_MAKEPED)
                    printf("\t%d", mother_id);
            }

            print_gender(person);
            print_name_id(person);
        }

        /* print extra individuals in tree for LINKAGE format */
        if (PRINTOUT_FORMAT == SHOW_PATHS)
        {
            for (temp = others_father; temp != NULL; temp = temp->next)
                printf("%s        N%d\t0\t0\t\t1\n", PEDIGREE_ID, temp->pid);
            for (temp = others_mother; temp != NULL; temp = temp->next)
                printf("%s        N%d\t0\t0\t\t2\n", PEDIGREE_ID, temp->pid);
        }

        if (PRINTOUT_FORMAT == FOR_MAKEPED)
        {
            for (temp = others_father; temp != NULL; temp = temp->next)
                if (member(temp->pid, individuals) == FALSE)
                    printf("%s        %d\t0\t0\t\t1\n", PEDIGREE_ID, temp->pid);
            for (temp = others_mother; temp != NULL; temp = temp->next)
                if (member(temp->pid, individuals) == FALSE)
                    printf("%s        %d\t0\t0\t\t2\n", PEDIGREE_ID, temp->pid);
        }
    }
}

/*In the graph described by root add in parent information to
  the records in the 2-d array info, in each row field 1 is
  father and field 2 is mother; count is number of people*/
void add_child_parent_link(info,root,count)
int **info;
tree_node *root;
int count;
{
    tree_node_list *temp; /*pointer into list of children*/
    int parent_gender;  /*gender of parent following LINKAGE conventions*/

    if ((root != NULL) && (root->children != NULL))
    {
        parent_gender = lookup_gender(root->self, root->children->this->self);
        /*Iterate over children*/
        for (temp = root->children; temp != NULL; temp = temp->next)
            add_link(root->self, temp->this->self, info, count, parent_gender);

        for (temp = root->children; temp != NULL; temp = temp->next)
           /*recursive call on child*/
            add_child_parent_link(info, temp->this, count);
    }
}

/*add information on parent_id into record for child_id,
  info is 2-d array of records, count is number of people, parent_gender
  is gender of parent_id*/
void add_link(parent_id, child_id, info, count, parent_gender)
int parent_id;
int child_id;
int **info;
int count;
int parent_gender;
{
    int index;  /*loop index over people */

    index = 0;

    while ((index < count) && (info[index][SELF] != child_id))
        index++;

    if (index == count)
    {
        fprintf(stderr, "Cannot find child %d in info\n", child_id);
        exit(EXIT_FAILURE);
    }

    if (parent_gender == MALE)
        info[index][FATHER_INDEX] = parent_id;
    else
        info[index][MOTHER_INDEX] = parent_id;
}

/*lookup the gender of person who has a child indexed by child*/
int lookup_gender(person, child)
int person;
int child;
{
    if (parent_information[child][FATHER] == person)
        return(MALE);

    if (parent_information[child][MOTHER] == person)
        return(FEMALE);

    fprintf(stderr,"%d not a parent of %d\n", person, child);
    exit(EXIT_FAILURE);
}


/*Return a list of individuals in the graph rooted at root*/
pid_list *list_individuals_down(root)
tree_node *root;
{
    tree_node_list *temp;  /*index over children of root*/
    pid_list *result = NULL; /*list to hold results*/

    if (root != NULL)
    {
        /*add self*/
        add_pid(root->self, &result);

        for (temp = root->children; temp != NULL; temp = temp->next)
          /*add list for each child*/
            union_list(&result, list_individuals_down(temp->this));
    }
    return(result);
}
