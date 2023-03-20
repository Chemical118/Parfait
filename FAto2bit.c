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

bool isfafile(char *f); //cheak the file is fasta
char* read_line(FILE *fp, char st); //read one line of the file and return it as a string
void read_fafile(FILE *fp, FILE *w_fp); //read the fasta file
int chr_to_num(char tmp); // ATGCN to 01234
char num_to_chr(int tmp); //01234 to ATGAN
Element mk_huffman_tree(Data_info data_tmp, char odd_data_tmp); //making huffman tree
void heap_insert(Heap *h, Element node); //insert node to heap
Element heap_delete(Heap *h); //delete node from heap
int encoding(Element head, Data_info data_tmp, FILE *w_fp); //encoding file with huffman tree
void destroy_heap(Heap_node *node);
void make_huffman_code(Heap_node *node, int len, char *code, char **code_arr);
void write_index_file(Data_info data_tmp);

int main(int argc, char *argv[]) {

    if (argc != 2) {
        printf("We need 1 FASTA File\n");
        return 0;
    }

    FILE *fasta_fp = NULL;
    fasta_fp = fopen(argv[1], "r");
    if (fasta_fp == NULL) {
        printf("File Reading ERROR\n");
        return 0;
    }

    if (isfafile(argv[1]) == false) {
        printf("Not a FASTA File\n");
        return 0;
    }

    FILE *print_fasta_fp = NULL;
//    print_fasta_fp = fopen(strcat(argv[1], ".fz"), "bw");
    print_fasta_fp = fopen("test1.fz", "wb");

    read_fafile(fasta_fp, print_fasta_fp);
    fclose(fasta_fp);
}

bool isfafile(char *f){
    char *point = NULL;
    point = strrchr(f, '.'); //last '.' in the file name
    if (strcmp(point, ".fasta") == 0) {
        return true;
    }
    return false;
}

char* read_line(FILE *fp, char st) {
    const int BUFF_SIZE = 100;
    int ch_cnt = 1;
    int buff_times = 1;
    char ch_tmp;
    char *line = (char*)malloc(BUFF_SIZE * buff_times * sizeof(char));
    line[0] = st;
    while( (ch_tmp = getc(fp)) != '\n') {
        line[ch_cnt] = ch_tmp;
        ch_cnt += 1;
        if (ch_cnt % 100 == 98) { //prevent overflow when adding '\n' and '\0'
            buff_times += 1;
            line = (char*)realloc(line, BUFF_SIZE * buff_times);
        }
    }
    line[ch_cnt++] = '\n';
    line[ch_cnt] = '\0';
    puts(line);
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

void heap_insert(Heap *h, Element node) {
    int index = h->size + 1;
    while ( (index != 1) && (node.freq < h->heap[index/2].freq) ) {
        h->heap[index] = h->heap[index/2];
        index = index/2;
    }
    h->heap[index] = node;
    h->size += 1;
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
        printf("%d %c%c : %d\n",idx,node->ch1,node->ch2,ele.freq);
    }

    //making huffman tree
    printf("heap size : %d\n",heap.size);
    int size = heap.size;
    for (int i = 0; i < size-1; i++) { //repeat until one node left

        //two minimum node
        e1 = heap_delete(&heap);
        e2 = heap_delete(&heap);
        printf("%d %d\n",e1.freq,e2.freq);
        //new node
        Heap_node *tmp = (Heap_node*)malloc(sizeof(Heap_node));
        tmp->left_node = e1.pnode;
        tmp->right_node = e2.pnode;
        Element ele;
        ele.freq = e1.freq + e2.freq;
        ele.pnode = tmp;

        heap_insert(&heap, ele);
    }

    return heap_delete(&heap); //return the head node of huffman tree
}

void make_huffman_code(Heap_node *node, int len, char *code, char **code_arr) {
    if (node != NULL) {

        len += 1;
        code[len] = '1';
        make_huffman_code(node->left_node, len, code, code_arr);

        code[len] = '0';
        make_huffman_code(node->right_node, len, code, code_arr);

        code[len] = '\0';

        if (node->left_node == NULL && node->right_node == NULL) {
            printf("%c%c\n",node->ch1,node->ch2);
            strcpy(code_arr[chr_to_num(node->ch1) * 5 + chr_to_num(node->ch2)], code);
            printf("code %s\n",code);
        }
    }
}

int encoding(Element head, Data_info data_tmp, FILE *w_fp) {
    int len = 0; //number of bits
    unsigned char buffer = '\0'; //8bit -> 1byte
    int bit_num = 0; //how many bits did we used at buffer
    if(w_fp == NULL) printf("!$!\n");

    char **code_arr = (char**)calloc(MAX_NODE, sizeof(char*)); //array that has huffman code
    for (int i = 0; i < MAX_NODE; i++) {
        code_arr[i] = (char*)calloc(100, sizeof(char));
    }
    char *code_str = (char*)calloc(100, sizeof(char));
    make_huffman_code(head.pnode, -1, code_str, code_arr);

    write_index_file(data_tmp);

    for(int i = 0; i < MAX_NODE; i++) {
        printf("%d : %s\n",i,code_arr[i]);
    }

    char *str = (char*)calloc(100, sizeof(char));
    for (int idx = 0; data_tmp.file_data[idx]; idx++) {
        strcpy(str, code_arr[data_tmp.file_data[idx]]);
        printf("%s", str);
        len += (int)strlen(str);
        for (int i = 0; str[i]; i++) {
            buffer = buffer << 1;
            buffer = buffer | (unsigned char) (str[i] - '0');
            bit_num += 1;

            if (bit_num == 8) {
                /*
                for(int j=0;j<8;j++) {
                    printf("%d", (buffer & (1 << j)) >> j);
                }
                printf("\n");
                */
                fwrite(&buffer, 1, 1, w_fp);
                buffer = '\0';
                bit_num = 0;

            }
        }
    }
    if (bit_num != 0) {
        while (bit_num < 8) {
            buffer = buffer << 1;
            bit_num += 1;
        }
        fwrite(&buffer, 1, 1, w_fp);
        /*
        for(int j=0;j<8;j++) {
            printf("%d", (buffer & (1 << j)) >> j);
        }
        printf("\n");
        */
    }

    printf("\n");

    free(code_arr);
    free(code_str);

    return len;
}

void read_fafile(FILE *fp, FILE *w_fp) {

    char ch_tmp;
    //data info

    Data_info data_tmp; //file data we read
    data_tmp.total_freq = 0;
    data_tmp.data_len = 0;
    data_tmp.data_num = 0;
    data_tmp.freq_arr = (int*)calloc(MAX_NODE, sizeof(int));
    data_tmp.file_data = (int*)calloc((MAX_LENGTH+10)/2, sizeof(int));

    int pre_data_num = 0;
    char odd_data_tmp = '\0'; //saving the first char

    while ((ch_tmp = getc(fp)) != EOF) { //all file reading
//        ch_tmp = getc(fp);
        if (ch_tmp == ';' || ch_tmp == '>') { //start of a new data
            data_tmp.data_name = read_line(fp, ch_tmp);
            if(data_tmp.data_len != 0) { //if we have some read data
                Element head = mk_huffman_tree(data_tmp, odd_data_tmp); //making huffman tree
                printf("%d end mk huffmantree\n",head.freq);
                encoding(head, data_tmp ,w_fp); //encoding the file
                printf("%d end encoding\n",data_tmp.data_len);
                destroy_heap(head.pnode);
                printf("end destroy heap\n");
                free(data_tmp.freq_arr);
                free(data_tmp.file_data);
                free(data_tmp.data_name);
                printf("end free\n");
                write_index_file(data_tmp);

                //initialization
                data_tmp.data_num = pre_data_num + 1;
//                data_tmp.data_name = read_line(fp, ch_tmp);
                data_tmp.freq_arr = (int*)calloc(MAX_NODE, sizeof(int));
                data_tmp.file_data = (int*)calloc((MAX_LENGTH+10)/2, sizeof(int));
                data_tmp.total_freq = 0;
                data_tmp.data_len = 0;
                pre_data_num += 1;
                printf("end init\n\n\n");
            }

            if (ch_tmp == ';') { //the line after symbol ';' doesn't need
                char *line_tmp = read_line(fp, ch_tmp);
                data_tmp.data_name = realloc(data_tmp.data_name,
                                             strlen(data_tmp.data_name) + strlen(line_tmp) + 1);
                strcat(data_tmp.data_name, line_tmp);
                free(line_tmp);
            }

            data_tmp.data_name[strlen(data_tmp.data_name) - 1] = '\0';


        } else { //reading the file

            if (ch_tmp == '\n' || ch_tmp == '\r' || ch_tmp == EOF) {

            } //don't read enter
            else {
                data_tmp.data_len += 1;
                if (data_tmp.data_len % 2 == 1) { //if it is the first char
                    odd_data_tmp = ch_tmp;
                } else { //if it is the second char
                    int index = chr_to_num(odd_data_tmp)*5+chr_to_num(ch_tmp);
                    odd_data_tmp = '\0';
                    data_tmp.freq_arr[index] += 1;
                    data_tmp.file_data[data_tmp.total_freq] = index; //recoding the file info
                    printf("file_data, idx:%d, %d\n",data_tmp.total_freq,index);
                    data_tmp.total_freq += 1;
                }
            }
        }

        if (data_tmp.data_len == MAX_LENGTH) { //if we read 'MAX_LENGTH' number of char
            Element head = mk_huffman_tree(data_tmp, odd_data_tmp); //make huffman tree
            encoding(head, data_tmp, w_fp); //encoding the file
            printf("%d end encoding",data_tmp.data_len);
            destroy_heap(head.pnode);
            free(data_tmp.freq_arr);
            free(data_tmp.file_data);
            free(data_tmp.data_name);
            write_index_file(data_tmp);

            //initialization
            data_tmp.freq_arr = (int*)calloc(MAX_NODE, sizeof(int));
            data_tmp.file_data = (int*)calloc((MAX_LENGTH+10)/2, sizeof(int));
            data_tmp.total_freq = 0;
            data_tmp.data_len = 0;
        }
    }

    if ( data_tmp.data_len != 0 ) {
        Element head = mk_huffman_tree(data_tmp, odd_data_tmp); //making huffman tree
        printf("%d end mk huffmantree\n",head.freq);
        encoding(head, data_tmp, w_fp); //encoding the file
        printf("%d end encoding\n",data_tmp.data_len);
        destroy_heap(head.pnode);
        printf("end destroy heap\n");
        free(data_tmp.freq_arr);
        free(data_tmp.file_data);
        free(data_tmp.data_name);
        printf("end free\n");
        write_index_file(data_tmp);
    }
}

void write_index_file(Data_info data_tmp) {
    printf("writing index file\n");
}
