#include "search.h"

char *get_token(char *, int);

void get_parent_information()
{
    int      pid1, pid2;
    char     *children;
    char     descriptions[MAXLINELEN];
    pid_list *children_pids = NULL;
    pid_list *temp;
    int      index;
    FILE     *infile;
    char     *token;

    infile = fopen(CHILD_PARENT_TABLE_NAME, "r"); 
    if (infile == NULL)
    { 
       fprintf(stderr, "Cannot open relationship table %s\n", 
                                        CHILD_PARENT_TABLE_NAME); 
       exit(EXIT_FAILURE);
    }

    parent_information = (int**) malloc((POPULATION+1) * sizeof(int *) );
    throw_if_null(parent_information,"parent_information");

    for (index = 0; index <= POPULATION; index++)
    {
        parent_information[index] = (int *) malloc(2 * sizeof(int));
        throw_if_null(parent_information[index], "parent_information");
    }

    couple_information = (int**) malloc((POPULATION+1) * sizeof(int *) );
    throw_if_null(couple_information, "couple_information");
    
    for (index = 0; index <= POPULATION; index++)
    {
        couple_information[index] = (int *) malloc(2 * sizeof(int));
        throw_if_null(couple_information[index], "couple_information");
    }

    for (index = 0; index <= POPULATION; index++)
    {
        parent_information[index][FATHER] = NO_RESULT;
        parent_information[index][MOTHER] = NO_RESULT;
        couple_information[index][FATHER] = NO_RESULT;
        couple_information[index][MOTHER] = NO_RESULT;
    }

    index = 1;
    while (fgets(descriptions, MAXLINELEN , infile) != NULL) 
    {
        /* Get father */
        token = get_token(descriptions, FATHER_COLUMN);
        if (token == NULL)
        {
            sprintf(error_message,"No father on line %d of relationship table\n",
                     index);
            print_exit(error_message);
        }
        pid1 = atoi(token);
        if ((pid1 < 1) || (pid1 > POPULATION))
        {
            sprintf(error_message,"%d is not a valid program id as it "
                    "is not in range 1..%d\n", pid1, POPULATION);
            print_exit(error_message);
        }

        /* Get mother */
        token = get_token(descriptions, MOTHER_COLUMN);
        if (token == NULL)
        {
            sprintf(error_message,"No mother on line %d of relationship table\n",
                     index);
            print_exit(error_message);
        }
        pid2 = atoi(token);
        if ((pid2 < 1) || (pid2 > POPULATION))
        {
            sprintf(error_message,"%d is not a valid program id as it "
                    "is not in range 1..%d\n", pid2, POPULATION);
            print_exit(error_message);
        }

        /* Get list of children */
        children = get_token(descriptions, CHILDREN_COLUMN);
        children_pids = NULL;
        list_for_children_ids(children, &children_pids);
        for (temp = children_pids; temp != NULL; temp = temp->next)
        {
            if ((parent_information[temp->pid][FATHER] != NO_RESULT) ||
                (parent_information[temp->pid][MOTHER] != NO_RESULT))
            {
                sprintf(error_message,"Error: %d present in more than one family\n", temp->pid);
                print_exit(error_message);
            }

            if ((temp->pid == pid1) || (temp->pid == pid2))
            {
                sprintf(error_message,"Error: %d is a parent and a child in the same family\n", temp->pid);
                print_exit(error_message);
            }

            parent_information[temp->pid][FATHER] = pid1;
            parent_information[temp->pid][MOTHER] = pid2;
        }
        couple_information[index][FATHER] = pid1;
        couple_information[index][MOTHER] = pid2;
        index++;
        if (children_pids != NULL)
            free_list(children_pids);
    }
    PARENT_TABLE = TRUE;
}

int lookup_father(person)
int person;
{
    if (PARENT_TABLE)
        return(parent_information[person][FATHER]);
    else
    {
        sprintf(error_message,"Download parent table first\n");
        print_exit(error_message);
    }
}

int lookup_mother(person)
int person;
{
    if (PARENT_TABLE)
        return(parent_information[person][MOTHER]);
    else
    {
        sprintf(error_message,"Download parent table first\n");
        print_exit(error_message);
    }
}

pid_list *lookup_spouses_with_children(person)
int person;
{
    pid_list *all_spouses = NULL;
    int index;

    if (PARENT_TABLE)
    {
        for (index = 1; index <= POPULATION; index++)
        {
            if (parent_information[index][FATHER] == person) 
                add_pid(parent_information[index][MOTHER], &all_spouses);
            if (parent_information[index][MOTHER] == person) 
                add_pid(parent_information[index][FATHER], &all_spouses);
        }
    }
    else
    {
        sprintf(error_message,"Download parent table first\n");
        print_exit(error_message);
    }

    return(all_spouses);
}

pid_list *lookup_spouses(person)
int person;
{
    pid_list *all_spouses = NULL;
    int index;

    if (PARENT_TABLE)
    {
        for (index = 1; index <= POPULATION; index++)
        {
            if (couple_information[index][FATHER] == person) 
                add_pid(couple_information[index][MOTHER], &all_spouses);
            if (couple_information[index][MOTHER] == person) 
                add_pid(couple_information[index][FATHER], &all_spouses);
        }
    }
    else
    {
        sprintf(error_message,"Download parent table first\n");
        print_exit(error_message);
    }

    return(all_spouses);
}

pid_list *lookup_children(person)
int person;
{
    pid_list *all_children = NULL;
    int index;

    if (PARENT_TABLE)
    {
        for (index = 1; index <= POPULATION; index++)
        {
            if ((parent_information[index][FATHER] == person) ||
                (parent_information[index][MOTHER] == person))
                add_pid(index, &all_children);
        }
    }
    else
    {
        sprintf(error_message,"Download parent table first\n");
        print_exit(error_message);
    }

    return(all_children);
}

pid_list *lookup_children_for_couple(person,spouse)
int person;
int spouse;
{
    pid_list *all_children = NULL;
    int index;

    if (PARENT_TABLE)
    {
        for (index = 1; index <= POPULATION; index++)
        {
            if (((parent_information[index][FATHER] == person) &&
                 (parent_information[index][MOTHER] == spouse)) ||
                ((parent_information[index][FATHER] == spouse)  &&
                 (parent_information[index][MOTHER] == person)))
                add_pid(index, &all_children);
        }
    }
    else
    {
        sprintf(error_message,"Download parent table first\n");
        print_exit(error_message);
    }

    return(all_children);
}

void get_person_information()
{
    int pid, index;
    FILE *infile;
    char *token;
    char     descriptions[MAXLINELEN];
    int lineNumber;

    infile = fopen(PERSON_TABLE_NAME, "r"); 
    if (infile == NULL)
    { 
       fprintf(stderr, "Cannot open person table %s\n", 
                                        PERSON_TABLE_NAME); 
       exit(EXIT_FAILURE);
    }

    person_information = (person_record *) malloc((POPULATION+1) * sizeof(person_record) );
    throw_if_null(person_information,"person_information");

    for (index = 0; index <= POPULATION; index++)
    {
        person_information[index].name = NULL;
        person_information[index].gender = 'U';  /* For unknown */
        person_information[index].birth_date = NULL;
        person_information[index].death_date = NULL;
        person_information[index].useful_information = NULL;
    }

    lineNumber = 1;
    while (fgets(descriptions, MAXLINELEN , infile) != NULL) 
    {
        /* Get pid */
        token = get_token(descriptions, PID_COLUMN);
        if (token == NULL)
        {
            sprintf(error_message,"No pid on line %d?\n", lineNumber);
            print_exit(error_message);
        }
        pid = atoi(token);
        if ((pid < 1) || (pid > POPULATION))
        {
            sprintf(error_message,"%d is not a valid program id as it "
                    "is not in range 1..%d\n", pid, POPULATION);
            print_exit(error_message);
        }

        /* Get name */
        person_information[pid].name = get_token(descriptions, NAME_COLUMN);
        
        /* Get birth_date */
        person_information[pid].birth_date = get_token(descriptions, BIRTH_DATE_COLUMN);

        /* Get death_date */
        person_information[pid].death_date = get_token(descriptions, DEATH_DATE_COLUMN);

        /* Get gender */
        token = get_token(descriptions, GENDER_COLUMN);
        if (token != NULL)
            person_information[pid].gender = token[0];
	person_information[pid].useful_information = get_token(descriptions, INFO_COLUMN);
	lineNumber++;
    }
    PERSON_TABLE = TRUE;
}

char *get_token(line,given_count)
char *line;
int  given_count;
{
    int index = 0;
    int count = 1;
    int len;
    char *result = NULL;
    int start_point = 0;
    int end_point = 0;

    if (line == NULL)
        return(NULL);
    len = strlen(line);

    /* find start point */
    while ((count < given_count) && (index < len))
    {
        if (line[index] == DELIMITER)
            count++;
        index++;
    }

    if (count == given_count)
        start_point = index;
    else
        return(NULL);

    /* find end point */
    while ((index < len) && (line[index] != DELIMITER))
        index++;
    end_point = index;

    if (end_point == start_point)
        return(NULL);

    result = (char *) malloc((end_point - start_point + 1) * sizeof(char));
    if (result == NULL)
    {
        sprintf(error_message,"Cannot allocate memory to get token %d\n",
                 given_count);
        print_exit(error_message);
    }

    for (index = start_point; index < end_point; index++)
        result[index-start_point] = line[index];
    if (end_point == len)
        result[index-start_point-1] = NIL;
    else
        result[index-start_point] = NIL;

    return(result);
}
