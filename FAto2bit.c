#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#define MAX_LENGTH 100000
#define MAX_NODE 25

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
    char ch1;
    char ch2;
}Heap_node; //A node from heap

typedef struct Element{
    int freq;
    Heap_node *pnode;
}Element; //data of heap node

typedef struct Heap{
    Element heap[MAX_NODE + 10];
    int size;
}Heap; //heap

bool isfafile(char *File); //cheak the file is fasta
char* read_line(FILE *fp); //read one line of the file and return it as a string
void read_fafile(FILE *fp); //read the fasta file
int chr_to_num(char tmp); // ATGCN to 01234
char num_to_chr(int tmp); //01234 to ATGAN
Element mk_huffman_tree(Data_info data_tmp, char odd_data_tmp); //making huffman tree
void heap_insert(Heap *h, Element node); //insert node to heap
Element heap_delete(Heap *h); //delete node from heap
void encoding(Element head, Data_info data_tmp); //encoding file with huffman tree
void destroy_heap(Heap_node *node);
void making_huffman_code(Heap_node *node, int len, char *code, int *code_arr);
int str_to_num(char *code);
char *num_to_str(int num);
void write_index_file(Data_info data_tmp);

int main(int argc, char *argv[]) {

    if (argc != 2) {
        printf("We need 1 FASTA File\n");
        return 0;
    }

    if (isfafile(argv[1]) == false) {
        printf("Not a FASTA File\n");
        return 0;
    }

    FILE *fasta_fp = NULL;
    fasta_fp = fopen(argv[1], "r");
    if (fasta_fp == NULL) {
        printf("File Reading ERROR\n");
        return 0;
    }
    read_fafile(fasta_fp);

}

bool isfafile(char *File) {
    if (strlen(File) % 2 == 0) return true;
    return false;
}

char* read_line(FILE *fp) {
    char ch_tmp;
    char *line = NULL;
    while( (ch_tmp = getc(fp)) != '\n') {
        char line_tmp[1] = {ch_tmp};
        strcat(line, line_tmp);
    }
    return line;
}

int chr_to_num(char tmp) {
    if (tmp == 'A' || tmp == 'a') {
        return 0;
    } else if (tmp == 'T' || tmp == 't') {
        return 1;
    } else if (tmp == 'G' || tmp == 'g') {
        return 2;
    } else if (tmp == 'C' || tmp == 'c') {
        return 3;
    } else if (tmp == 'N' || tmp == 'n'){
        return 4;
    } else return -1;
}

char num_to_chr(int tmp) {
    if (tmp == 0) {
        return 'A';
    } else if (tmp == 1) {
        return 'T';
    } else if (tmp == 2) {
        return 'G';
    } else if (tmp == 3) {
        return 'C';
    } else if (tmp == 4){
        return 'N';
    } else return '\0';
}

int str_to_num(char *code) {
    int num = 0;
    int k = 1;
    for (int i = (int)strlen(code) - 1; i >= 0; i--) {
        num += ((code[i] - '0') * k);
        k *= 2;
    }
    return num;
}

char *num_to_str(int num) {
    int len = 0;
    char *str = (char*)calloc(100, sizeof(char));
    for (int i = 1; i <= num; i*=2) { //calculate the length of the string
        len += 1;
    }
    for (int i = 0; i < len; i++) {
        if (((1 << (len - i)) & num) == 1) {
            str[i] = '1';
        } else {
            str[i] = '0';
        }
    }
    return str;
}

void heap_insert(Heap *h, Element node) {
    int index = h->size + 1;
    while ( (index != 1) && (node.freq < h->heap[index/2].freq) ) {
        h->heap[index] = h->heap[index/2];
        index = index/2;
    }
    h->heap[index] = node;
}

Element heap_delete(Heap *h) {
    int parent = 1, child = 2;
    Element node, tmp;

    node = h->heap[1];
    tmp = h->heap[h->size];
    h->size -= 1;

    while (child < h->size) {
        if (child != h->size && ( (h->heap[child].freq) > (h->heap[child+1].freq) )) {
            child += 1;
        }

        if (tmp.freq <= h->heap[child].freq) {
            break;
        }

        h->heap[parent] = h->heap[child];
        parent = child;
        child = child * 2;
    }

    h->heap[parent] = tmp;
    return node;
}

void destroy_heap(Heap_node *node) {
    if (node == NULL) return;

    destroy_heap(node->left_node);
    destroy_heap(node->right_node);
    free(node);
}

Element mk_huffman_tree(Data_info data_tmp, char odd_data_tmp) {
    Heap heap; //min heap
    heap.size = 0;
    Element e1,e2; //minimum 2 elements

    if (odd_data_tmp != '\0') { //if the data is odd
        int index = chr_to_num(odd_data_tmp) * 5 + chr_to_num('A'); //add 'A' at the end
        data_tmp.freq_arr[index] += 1;
        data_tmp.file_data[data_tmp.total_freq] = index;
        data_tmp.total_freq += 1;
    }

    for (int idx = 0; idx < MAX_NODE; idx++) { //first heap node
        if (data_tmp.freq_arr[idx] == 0) continue; //if it doesn't appear at file

        Heap_node *node = (Heap_node*)malloc(sizeof(Heap_node)); //new node
        Element ele; //node's info element
        ele.freq = data_tmp.freq_arr[idx];
        node->ch1 = num_to_chr(idx/5); //first char
        node->ch2 = num_to_chr(idx%5); //second char
        ele.pnode = node; //Element -> Heap node
        //it is leaf node
        node->left_node = NULL;
        node->right_node = NULL;

        heap_insert(&heap, ele);
    }

    //making huffman tree
    for (int i = 0; i < heap.size-1; i++) { //repeat until one node left

        //two minimum node
        e1 = heap_delete(&heap);
        e2 = heap_delete(&heap);

        //new node
        Heap_node *tmp = (Heap_node*)malloc(sizeof(Heap_node));
        tmp->right_node = e1.pnode;
        tmp->left_node = e2.pnode;
        Element ele;
        ele.freq = e1.freq + e2.freq;
        ele.pnode = tmp;

        heap_insert(&heap, ele);
    }

    return heap_delete(&heap); //return the head node of huffman tree
}

void making_huffman_code(Heap_node *node, int len, char *code, int *code_arr) {
    if (node != NULL) {
        len += 1;
        code[len] = '1';
        making_huffman_code(node->left_node,len,code,code_arr);

        code[len] = '0';
        making_huffman_code(node->right_node,len,code,code_arr);

        code[len] = '\0';

        if (node->left_node == NULL && node->right_node == NULL) {
            code_arr[chr_to_num(node->ch1) * 5 + chr_to_num(node->ch2)] = str_to_num(code);
        }
    }
}

void encoding(Element head, Data_info data_tmp) {
    int *code_arr = (int*)calloc(MAX_NODE, sizeof(int)); //array that has huffman code
    char *code_str = (char*)calloc(100, sizeof(char));
    making_huffman_code(head.pnode, -1, code_str, code_arr);

    write_index_file(data_tmp);

    for (int idx = 0; data_tmp.file_data[idx]; idx++) {
        printf("%s", num_to_str(code_arr[data_tmp.file_data[idx]]));
    }

    free(code_arr);
    free(code_str);
}

void read_fafile(FILE *fp) {

    char ch_tmp;
    //data info

    Data_info data_tmp; //file data we read
    int pre_data_num = 0;
    char odd_data_tmp = '\0'; //saving the first char


    while ( (ch_tmp = getc(fp)) != EOF) { //all file reading
        if (ch_tmp == ';' || ch_tmp == '>') { //start of a new data
            if(data_tmp.total_freq != 0) { //if we have some read data
                Element head = mk_huffman_tree(data_tmp, odd_data_tmp); //making huffman tree
                encoding(head, data_tmp); //encoding the file
                destroy_heap(head.pnode);
                free(data_tmp.freq_arr);
                free(data_tmp.file_data);

                //initialization
                data_tmp.data_num = pre_data_num + 1;
                data_tmp.data_name = read_line(fp);
                data_tmp.freq_arr = (int*)calloc(MAX_NODE, sizeof(int));
                data_tmp.file_data = (int*)calloc((MAX_LENGTH+10)/2, sizeof(int));
                data_tmp.total_freq = 0;
                data_tmp.data_len = 0;
                pre_data_num += 1;
            }

            if (ch_tmp == ';') { //the line after symbol ';' doesn't need
                read_line(fp);
            }

        } else { //reading the file
            if (ch_tmp == '\n') continue; //don't read enter
            else {
                data_tmp.data_len += 1;
                if (data_tmp.data_len % 2 == 1) { //if it is the first char
                    odd_data_tmp = ch_tmp;
                } else { //if it is the second char
                    int index = chr_to_num(odd_data_tmp)*5+chr_to_num(ch_tmp);
                    odd_data_tmp = '\0';
                    data_tmp.freq_arr[index] += 1;
                    data_tmp.file_data[data_tmp.total_freq] = index; //recoding the file info
                    data_tmp.total_freq += 1;
                }
            }
        }

        if (data_tmp.data_len == MAX_LENGTH) { //if we read 'MAX_LENGTH' number of char
            Element head = mk_huffman_tree(data_tmp, odd_data_tmp); //make huffman tree
            encoding(head, data_tmp); //encoding the file
            destroy_heap(head.pnode);
            free(data_tmp.freq_arr);
            free(data_tmp.file_data);

            //initialization
            data_tmp.freq_arr = (int*)calloc(MAX_NODE, sizeof(int));
            data_tmp.file_data = (int*)calloc((MAX_LENGTH+10)/2, sizeof(int));
            data_tmp.total_freq = 0;
            data_tmp.data_len = 0;
        }
    }
}