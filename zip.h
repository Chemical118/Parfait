#ifndef PARFAIT_ZIP_H
#define PARFAIT_ZIP_H

typedef struct Data_info{
    int data_num;
    char *data_name;
    int *freq_arr;
    int *file_data;
    int total_freq;
    int data_len;
}Data_info; //data about the fasta file

typedef struct Heap_node{
    struct Heap_node *left_node;
    struct Heap_node *right_node;
    char *ch;
}Heap_node; //A node from heap

typedef struct Element{
    int freq;
    Heap_node *pnode;
}Element; //data of heap node

typedef struct Heap{
    Element *heap;
    int size;
}Heap; //heap

#endif //PARFAIT_ZIP_H
