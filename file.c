#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *load_file_content(char *file_name)
{
    char *file_contents;
    long input_file_size;

    FILE *input_file = fopen(file_name, "r");
    if (input_file == NULL)
    {
        printf("Error in function file_reader.c/load_file_content: File not found\n");
        exit(1);
    }

    fseek(input_file, 0, SEEK_END);
    input_file_size = ftell(input_file);
    rewind(input_file);

    file_contents = malloc((input_file_size + 1) * (sizeof(char)));
    if(fread(file_contents, sizeof(char), input_file_size, input_file) != input_file_size){
        printf("Error in function file_reader.c/load_file_content: File not found\n");
        exit(1);
    };
    file_contents[input_file_size] = '\0';

    fclose(input_file);

    return file_contents;
}

char *get_key_value(char *root, char *key)
{
    char *root_copy;
    int key_len = strlen(key);
    char *begin_tag_key = (char*) malloc(sizeof(char) * (key_len + 2));
    char *end_tag_key = (char*) malloc(sizeof(char) * (key_len + 3));
    begin_tag_key[0]  = '<';
    begin_tag_key[1]  = '\0';
    end_tag_key[0] = '<';
    end_tag_key[1] = '/';
    end_tag_key[2] = '\0';

    root_copy = (char*) malloc(sizeof(char) * (strlen(root) + 1));
    strcpy(root_copy, root);

    char *begin_tag_pointer = strstr(root_copy, strcat(begin_tag_key, key));
    char *end_tag_pointer = strstr(root_copy, strcat(end_tag_key, key));
    if (begin_tag_pointer == NULL || end_tag_pointer == NULL)
    {
      printf("Error in function file_reader.c/load_file_content: Key %s not found in the file\n", key);
      exit(1);
    }

    end_tag_pointer[0] = '\0';
    int return_value_len = strlen(begin_tag_pointer + key_len + 2);
    char *return_value = (char*) malloc(sizeof(char) * (1 + return_value_len));
    strcpy(return_value, begin_tag_pointer + key_len + 2);
    if (strlen(return_value) < 1) {
        printf("Error in function file_reader.c/load_file_content: Invalid value for key %s\n", key);
        exit(1);
    }

    free(begin_tag_key);
    free(end_tag_key);
    free(root_copy);
    return return_value;
}

double char_to_double(char *c) {
  double return_value = strtod(c, NULL);
  free(c);
  return return_value;
}

double char_to_int(char *c) {
  double return_value = atoi(c);
  free(c);
  return return_value;
}
