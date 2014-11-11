/*
 * Process lists, keep them in sorted order
*/ 


#include "search.h"

/*
 * throw_if_null() used for debugging purposes
 */

void throw_if_null(void* ptr, char* name)
{
    if (ptr == NULL)
    {
        fprintf(stderr,"No space for %s\n", name);
        exit(EXIT_FAILURE);
    }
}

int length_list(given_list)
pid_list *given_list;
{
    pid_list *temp;
    int count = 0;

    for (temp = given_list; temp != NULL; temp = temp->next)
        count++;
    return(count);
}

void print_list(given_list)
pid_list *given_list;
{
    pid_list *current = given_list;
    int count = 0;

    while (current != NULL)
    {
        if (count == LENGTH_CUTOFF)
        {
           printf("\n");
           count = 0;
        }
        printf("%d ", current->pid);
        current = current->next;
        count++;
    }
    printf("\n");
}

void free_list(given_list)
pid_list *given_list;
{
    pid_list *current = given_list;
    pid_list *remainder;

    while (current != NULL)
    {
        remainder = current->next;
        free(current);
        current = remainder;
    }
}

/* return list1 - list2 */
pid_list *subtract_list(list1, list2)
pid_list *list1;
pid_list *list2;
{
      pid_list *temp;
      pid_list *temp2;
      pid_list *result;

      result = duplicate_list(list1);

      while ((result != NULL) && (member(result->pid, list2)))
      {
         temp = result->next;
         free(result);
         result = temp;
      }

      if (result == NULL)
         return(NULL);

      /* First pid in result is not in list2 */

      for (temp = result; temp->next != NULL; )
           if (member(temp->next->pid, list2))
           {
                temp2 = temp->next;
                temp->next = temp->next->next;
                free(temp2);
           }
           else
               temp = temp->next;

      return(result);
}

void delete_pid(given_list, pid_to_delete)
pid_list **given_list;
int pid_to_delete;
{
    pid_list *current = *given_list;
    pid_list *previous;

    if ((current != NULL) && (current->pid == pid_to_delete))
    {
       *given_list = current->next;
       free(current);
    }
    else
    {
       while ((current != NULL) && (current->pid != pid_to_delete))
       {  
           previous = current;
           current = current->next;
       }

       if (current != NULL)  /* We are at the node to be deleted */
       {
           previous->next = current->next;
           free(current);
       }
    }
}

void add_pid(given_pid, given_list)
int given_pid;
pid_list **given_list;
{
    pid_list *new;
    pid_list *previous;

    if (!(member(given_pid, *given_list)))
    {
        new = (pid_list *) malloc(sizeof(pid_list));
        if (new == NULL)
        {
            sprintf(error_message,
		"add_pid: Can't allocate space for adding %d\n", given_pid);
            print_exit(error_message);
        }

        new->pid = given_pid;
        new->next = NULL;

        if (*given_list == NULL)
            *given_list = new;
        else
        {
            previous = *given_list; 
            if (previous->pid > given_pid)
            {
                new->next = *given_list;
                *given_list = new;
            }
            else
            {
                while ((previous->next != NULL) && 
                       (previous->next->pid < given_pid))
                    previous = previous->next;

                new->next = previous->next;
                previous->next = new; 
            }
        }
    }
}

/* append pid to beginning of list */
void append_pid(given_pid, given_list)
int given_pid;
pid_list **given_list;
{
    pid_list *new;

    new = (pid_list *) malloc(sizeof(pid_list));
    if (new == NULL)
    {
        sprintf(error_message,
		"append_pid: Can't allocate space for adding %d\n", given_pid);
        print_exit(error_message);
    }

    new->pid = given_pid;
    new->next = *given_list;
    *given_list = new;
}

pid_list *duplicate_list(given_list)
pid_list *given_list;
{
    pid_list *result = NULL;
    pid_list *current;

    for (current = given_list; current != NULL; current = current->next)
        append_pid(current->pid, &result);

    /* list is in reverse sorted order */
    reverse_list(&result);

    return(result);
}

void reverse_list(given_list)
pid_list **given_list;
{
     pid_list *current;
     pid_list *temp;
     pid_list *result = NULL;

     for (current = *given_list; current != NULL; )
     {
        append_pid(current->pid, &result);
        temp = current;
        current = current->next;
        free(temp);
     }

     *given_list = result;
}

/* We know that given_list is sorted in ascending order */
bool member(given_pid, given_list)
int given_pid;
pid_list *given_list;
{
    pid_list *current = given_list;

    while ((current != NULL) && (current->pid < given_pid))
       current = current->next;

    if ((current != NULL) && (current->pid == given_pid))
       return(TRUE);
    else
       return(FALSE);
}

void union_list(result, list2)
pid_list **result;
pid_list *list2;
{
    pid_list *current;

    for (current = list2; current != NULL; current = current->next)
        add_pid(current->pid, result);
}

pid_list *intersection_list(list1, list2)
pid_list *list1;
pid_list *list2;
{
    pid_list *result = NULL;
    pid_list *current;

    for (current = list1; current != NULL; current = current->next)
       if(member(current->pid, list2))
          add_pid(current->pid, &result);

    return(result);
}

pid_list *read_pids(input_file)
char *input_file;
{
    FILE *infile;
    pid_list *input_pids = NULL;
    int  current;
    char nextline[MAXLINELEN];

    infile = fopen(input_file,"r");

    if (infile == NULL)
    {
        sprintf(error_message,"Can't open input file %s\n", input_file);
        print_exit(error_message);
    }

    while (fgets(nextline, MAXLINELEN, infile) != NULL)
    {
        sscanf(nextline,"%d",&current);
        add_pid(current, &input_pids);
    }

    fclose(infile);
    return(input_pids);
}

/* 
   Form a list of pids in formatted_pids from all_pids
   where all_pids is of the form 
      <> id_1 <><> id_2 <> ... <> id_n <>
*/
void list_for_children_ids(all_pids, formatted_pids)
char *all_pids;
pid_list **formatted_pids;
{
    char *token;

    if (all_pids != NULL)
    {
        /* Use the fact that all_pids is of form
              <> id_1 <><> id_2 <> ... <> id_n <> */

        token = (char *) strtok(all_pids, WHITE_SPACE); /* get <> */
        while (token != NULL)
        {
            token = (char *) strtok(NULL, WHITE_SPACE);     /* get id_x */
      
            if (token != NULL)
            {
                add_pid(atoi(token), formatted_pids);
                token = (char *) strtok(NULL, WHITE_SPACE); /* get <> or <><> */
            }
        }
    }
}
