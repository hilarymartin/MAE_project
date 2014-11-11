/*
 * dbaccess.c
 *
 * Routines for interfacing with database
 */

#include "search.h"

void parse_arguments(number_of_arguments, argument_values)
int   number_of_arguments;
char  **argument_values;
{
    int current_arg;

    /* defaults for variable names */
    PERSON = 0;
    PERSON2 = 0;
    INFILE_NAME = NULL;
    OUTFILE_NAME = NULL;
    PARENT_TABLE = FALSE;
    PERSON_TABLE = FALSE;

    for (current_arg = 1; current_arg < number_of_arguments; current_arg++)
    {
        if (argument_values[current_arg][0] != '-')
        {
            printf("Invalid argument: %s\n", argument_values[current_arg]);
            exit(EXIT_FAILURE);
        }

        switch(argument_values[current_arg][1]) 
        {
             case 'X':  int_suffix(&PERSON, argument_values[current_arg],2);
                        break;
             case 'Y':  int_suffix(&PERSON2, argument_values[current_arg],2);
                        break;
             case 'F':  suffix(&INFILE_NAME, argument_values[current_arg], 2);
                        break;
             case 'K':  int_suffix(&KEY, argument_values[current_arg], 2);
                        break;
             case 'O':  suffix(&OUTFILE_NAME, argument_values[current_arg], 2);
                        break;
             default:   printf("Invalid argument: %s -- Ignoring it\n", 
                                       argument_values[current_arg]);
        }
    }
}

/* Store suffix of argument_value after offset in variable_name */
void suffix(variable_name, argument_value, offset)
char **variable_name;
char *argument_value;
int offset;
{
    int len = strlen(argument_value);
    char *result;

    if (len == offset)
    {
       *variable_name = NULL;
       return;
    }

    result = (char *) malloc((len - offset + 1) * sizeof(char));
    if (result == NULL)
    {
       printf("Can't allocate space for storing command line arguments\n");
       exit(EXIT_FAILURE);
    }

    result[len - offset] = NIL;
    for ( ; len >= offset; len--)
       result[len - offset - 1] = argument_value[len - 1];

    *variable_name = result;
}

/* Store int suffix of argument_value after offset in variable_name */
void int_suffix(variable_name, argument_value, offset)
int  *variable_name;
char *argument_value;
int offset;
{
   char *result;

   suffix(&result, argument_value, offset);
   if (result == NULL)
    {
        fprintf(stderr,"No integer value after %s\n", argument_value);
        exit(EXIT_FAILURE);
    }
   *variable_name = atoi(result);
}

bool correct_argument(flag)
int flag;
{
    switch(flag) {
       case 'X': if (PERSON == 0)
                     return(FALSE); 
                 if (PERSON > POPULATION)
                 {
                    fprintf(stderr,"Program id of a person cannot be greater than %d\n", POPULATION);
                    return(FALSE); 
                 }
                 if (PERSON < 0)
                 {
                    fprintf(stderr,"Program id of a person cannot be negative\n");
                    return(FALSE); 
                 }
                 break;
       case 'Y': if (PERSON2 == 0)
                     return(FALSE); 
                 if (PERSON2 > POPULATION)
                 {
                     fprintf(stderr,"Program id of a person cannot be greater than %d\n", POPULATION);
                     return(FALSE); 
                 }
                 if (PERSON2 < 0)
                 {
                    fprintf(stderr,"Program id of a person cannot be negative\n");
                    return(FALSE); 
                 }
                 break;
       case 'K': if ((KEY < 0) || (KEY > 2)) 
		 {
		    fprintf(stderr,"Key must be 0(all), 1(males only), or 2(females only)");
		    return(FALSE);
		 }
		 break;
       case 'F': if (INFILE_NAME == NULL)
                     return(FALSE); 
                 break;
       case 'O': if (OUTFILE_NAME == NULL)
                     return(FALSE); 
                 break;
       default: return(FALSE);
    }

    return(TRUE);
}

void print_exit(message)
char *message;
{
    fprintf(stderr, message);
    exit(EXIT_FAILURE);
}
