#include "search.h"

static int gender(int);

void print_gender(person)
int person;
{
     switch(gender(person)) {
         case MALE:    printf("\t\t1"); break;
         case FEMALE:  printf("\t\t2"); break;
         default:      printf("\t\t0"); 
     }
}

void print_name_id(person)
int person;
{
    char *description, *first_name;

    description = get_name(person);
    if (description == NULL)
        printf("\t// <name unknown> ");
    else
    {
        first_name = strtok(description, WHITE_SPACE);
        if (first_name == NULL)
            printf("\t// <name unknown> ");
        else
            printf("\t// %s ", first_name);
    }
 
    description = get_id(person);
    if (description != NULL)
        printf("(%s)\n", description);
    else
        printf("\n");
}

char *get_id(person)
int person;
{
    return(NULL);
}

int gender(person)
int person;
{
    if (PERSON_TABLE)
    {
        if ((person < 0) || (person > POPULATION))
            return(UNKNOWN_GENDER);
        if (person_information[person].gender == 'M')
            return(MALE);
        if (person_information[person].gender == 'F')
            return(FEMALE);
        return(UNKNOWN_GENDER);
    }
    else
    {
        sprintf(error_message,"Download person table first\n");
        print_exit(error_message);
    }
}

char *get_name(person)
int person;
{
    if (PERSON_TABLE)
    {
        if ((person < 0) || (person > POPULATION))
            return(NULL);
        else
            return(person_information[person].name);
    }
    else
    {
        sprintf(error_message,"Download person table first\n");
        print_exit(error_message);
    }
}

void print_person_information(person)
int person;
{
    char *description;

    if (!(PERSON_TABLE))
    {
        sprintf(error_message,"Download person table first\n");
        print_exit(error_message);
    }

    if (person_information[person].name == NULL)
        printf("<name unknown>");
    else
        printf("%s", person_information[person].name);
 
    description = get_id(person);
    if (description != NULL)
        printf(" (%s)", description);

    if (person_information[person].birth_date != NULL)
        printf("\t%s", person_information[person].birth_date);
    else
        printf("\t00/00/0000");

    if (person_information[person].death_date != NULL)
        printf("\t%s", person_information[person].death_date);
    else
        printf("\t00/00/0000");

    printf("\n");
}

void print_rel_linkage(individuals, visited)
pid_list *individuals;
char *visited;
{
    pid_list *temp;
    int person, father_id, mother_id;

    if (!(PERSON_TABLE))
    {
        sprintf(error_message,"Download person table first\n");
        print_exit(error_message);
    }

    for (temp = individuals; temp != NULL; temp = temp->next)
    {
        person = temp->pid;
        printf("%s\tP%d", PEDIGREE_ID, person);

        father_id = lookup_father(person); 
        if (father_id <= NO_RESULT)
            printf("\t0");
        else
        {
          if (visited[father_id])
	    printf("\tP%d", father_id);
          else
          {
	    printf("\tN%d", father_id);
            visited[father_id] = 2;
          }
        }

        mother_id = lookup_mother(person); 
        if (mother_id <= NO_RESULT)
            printf("\t0");
        else
        {
          if (visited[mother_id])
	    printf("\tP%d", mother_id);
          else
          {
	    printf("\tN%d", mother_id);
            visited[mother_id] = 2;
          }
        }

        print_gender(person);
        print_name_id(person);
    }

    for (person=0; person<=POPULATION; person++)
        if (visited[person] == 2)
        {
            printf("%s\tN%d\t0\t0", PEDIGREE_ID, person);
            print_gender(person);
            printf("\n");
        }
}

void print_list_linkage(individuals)
pid_list *individuals;
{
    pid_list *temp;
    pid_list *others;
    int person, father_id, mother_id;

    /* parents which are not in result graph but have to be in pedigree
       output because of LINKAGE format */
    others = NULL;

    for (temp = individuals; temp != NULL; temp = temp->next)
    {
        person = temp->pid;
        printf("%s\tP%d", PEDIGREE_ID, person);

        father_id = lookup_father(person); 
        if (father_id <= NO_RESULT)
            printf("\t0");
        else
        {
            if (member(father_id, individuals))
                printf("\tP%d", father_id);
            else
            {
                add_pid(father_id, &others);
                printf("\tN%d", father_id);
            }
        }

        mother_id = lookup_mother(person); 
        if (mother_id <= NO_RESULT)
            printf("\t0");
        else
        {
            if (member(mother_id, individuals))
                printf("\tP%d", mother_id);
            else
            {
                add_pid(mother_id, &others);
                printf("\tN%d", mother_id);
            }
        }

        print_gender(person);
        print_name_id(person);
    }

    /* print extra individuals in tree for LINKAGE format */
    for (temp = others; temp != NULL; temp = temp->next)
    {
        printf("%s\tN%d\t0\t0", PEDIGREE_ID, temp->pid);
        print_gender(temp->pid);
        printf("\n");
    }
}
