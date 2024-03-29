#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include "zip.h"

const int MAX_LENGTH = 100000;
//const int MAX_NODE = 25;

static int power5(int power);
static bool isfafile(char *f); //cheak the file is fasta
static char* read_line(FILE *fp, char st); //read one line of the file and return it as a string
static void read_fafile(FILE *fp, FILE *w_fp, int tmp_len); //read the fasta file
static int chr_to_num(char tmp);
static int str_to_num(const char *str, int tmp_len); // ATGCN to 01234
static char* num_to_chr(int idx, int tmp_len); //01234 to ATGAN
static Element mk_huffman_tree(Data_info data_tmp, char *odd_data_tmp, int tmp_len); //making huffman tree
static void heap_insert(Heap *h, Element node); //insert node to heap
static Element heap_delete(Heap *h); //delete node from heap
static int encoding(Element head, Data_info data_tmp, FILE *w_fp, int tmp_len); //encoding file with huffman tree
static void destroy_heap(Heap_node *node);
static void make_huffman_code(Heap_node *node, int len, char *code, char **code_arr, int tmp_len);
static void write_index_file(Data_info data_tmp);
int main_zip(int argc, char *argv[]);

// verbose level
static int verbose_level = 0;

static int usage(FILE *fp, int status) {
    fprintf(fp, "test\n");
    return status;
}


int main_zip(int argc, char *argv[]) {
    static const struct option loptions[] = {
        {"help", no_argument, NULL, 'h'},
        {"force", no_argument, NULL, 'f'},
        {"reindex", no_argument, NULL, 'r'},
        {"keep", no_argument, NULL, 'k'},
        {"verbose", no_argument, NULL, 'v'},
        {"encode_length", required_argument, NULL, 'l'},
        {NULL, no_argument, NULL, 0}
    };

    int c, encode_length = 3;
    bool is_force = false, reindex = false, keep = false;

    if (argc == 1) return usage(stderr, 1);

    while ((c = getopt_long(argc, argv, "hfrkvl:", loptions, NULL)) >= 0) {
        switch (c) {
            case 'h': return usage(stdout, 0);
            case '?': return usage(stdout, 1);
            case 'l': encode_length = atoi(optarg); break;
            case 'f': is_force = true; break;
            case 'k': keep = true; break;
            case 'v': verbose_level += 1; break;
        }
    }

    FILE *fasta_fp = NULL;
    fasta_fp = fopen(argv[optind], "r");
    if (fasta_fp == NULL) {
        fprintf(stderr, "Could not open file: %s\n", argv[optind]);
        return 1;
    }

    FILE *print_fasta_fp = NULL;
    char *filename = (char*)malloc((strlen(argv[optind]) + 5) * sizeof(char));
    strcpy(filename, argv[optind]);
    strcat(filename, ".fz");
    print_fasta_fp = fopen(filename, "wb");

    read_fafile(fasta_fp, print_fasta_fp, encode_length);
    fclose(fasta_fp);
    return 0;
}

static bool isfafile(char *f) {
    char *point = NULL;
    point = strrchr(f, '.'); //last '.' in the file name
    if (strcmp(point, ".fasta") == 0 || strcmp(point, ".fa") == 0) {
        return true;
    }
    return false;
}

static char* read_line(FILE *fp, char st) {
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
    //puts(line);
    return line;
}

static int power5(int power) {
    int value = 1;
    for (int i = 0; i < power; i++) {
        value *= 5;
    }
    return value;
}

static int chr_to_num(char tmp) {
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

static int str_to_num(const char *str, int tmp_len) {
    int value = 0;
    int k = 1;
    for (int i = tmp_len - 1; i >= 0; i--) {
        char tmp = str[i];
        value += chr_to_num(tmp) * k;
        if (chr_to_num(tmp) == -1) {
            return -1;
        }
        k *= 5;
    }
    return value;
}

static char* num_to_chr(int idx, int tmp_len) {
    char *value = (char*)calloc((tmp_len + 2), sizeof(char));
    int k = 1;
    for (int i = 1; i < tmp_len; i++, k *= 5);
    for (int i = 1; i <= tmp_len; i++, k /= 5) {
        int num = idx/k;
        idx = idx - num * k;
        if (num == 0) {
            value[i] = 'A';
        } else if (num == 1) {
            value[i] = 'T';
        } else if (num == 2) {
            value[i] = 'G';
        } else if (num == 3) {
            value[i] = 'C';
        } else if (num == 4) {
            value[i] = 'N';
        } else {
            return NULL;
        }
    }
    return value + 1;
}

static void heap_insert(Heap *h, Element node) {
    int index = h->size + 1;
    while ( (index != 1) && (node.freq < h->heap[index/2].freq) ) {
        h->heap[index] = h->heap[index/2];
        index = index/2;
    }
    h->heap[index] = node;
    h->size += 1;
}

static Element heap_delete(Heap *h) {
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

static void destroy_heap(Heap_node *node) {
    if (node == NULL) return;

    destroy_heap(node->left_node);
    destroy_heap(node->right_node);
    free(node);
}

static Element mk_huffman_tree(Data_info data_tmp, char *odd_data_tmp, int tmp_len) {
    const int MAX_NODE = power5(tmp_len);

    Heap heap; //min heap
    heap.size = 0;
    heap.heap = (Element*)malloc((MAX_NODE + 10) * sizeof(Element));
    Element e1,e2; //minimum 2 elements
    //printf("##%d\n",(int) strlen(odd_data_tmp));
    if (strlen(odd_data_tmp) != 0) { //if the data is odd
        for (int i = (int)strlen(odd_data_tmp); i <= tmp_len; i++) {
            odd_data_tmp[i] = 'A';
        }
        int index = str_to_num(odd_data_tmp + 1, tmp_len); //add 'A' at the end
        data_tmp.freq_arr[index] += 1;
        data_tmp.file_data[data_tmp.total_freq] = index;
        data_tmp.total_freq += 1;
    }
    //printf("%s\n",odd_data_tmp + 1);

    for (int idx = 0; idx < MAX_NODE; idx++) { //first heap node
        if (data_tmp.freq_arr[idx] == 0) continue; //if it doesn't appear at file

        Heap_node *node = (Heap_node*)malloc(sizeof(Heap_node)); //new node
        node->ch = (char*)calloc((tmp_len + 1), sizeof(char));
        Element ele; //node's info element
        ele.freq = data_tmp.freq_arr[idx];
        node->ch = num_to_chr(idx, tmp_len);
        ele.pnode = node; //Element -> Heap node
        //it is leaf node
        node->left_node = NULL;
        node->right_node = NULL;

        heap_insert(&heap, ele);
        //printf("%d %s : %d\n",idx,node->ch,ele.freq);
    }

    //making huffman tree
    //printf("heap size : %d\n",heap.size);
    int size = heap.size;
    for (int i = 0; i < size-1; i++) { //repeat until one node left

        //two minimum node
        e1 = heap_delete(&heap);
        e2 = heap_delete(&heap);
        //printf("%d %d\n",e1.freq,e2.freq);
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

static void make_huffman_code(Heap_node *node, int len, char *code, char **code_arr, int tmp_len) {
    if (node != NULL) {

        len += 1;
        code[len] = '1';
        make_huffman_code(node->left_node, len, code, code_arr, tmp_len);

        code[len] = '0';
        make_huffman_code(node->right_node, len, code, code_arr, tmp_len);

        code[len] = '\0';

        if (node->left_node == NULL && node->right_node == NULL) {
            //printf("%s\n",node->ch);
            strcpy(code_arr[str_to_num(node->ch, tmp_len)], code);
            //printf("code %s\n",code);
        }
    }
}

static int encoding(Element head, Data_info data_tmp, FILE *w_fp, int tmp_len) {
    const int MAX_NODE = power5(tmp_len);
    int len = 0; //number of bits
    unsigned char buffer = '\0'; //8bit -> 1byte
    int bit_num = 0; //how many bits did we used at buffer

    char **code_arr = (char**)calloc(MAX_NODE, sizeof(char*)); //array that has huffman code
    for (int i = 0; i < MAX_NODE; i++) {
        code_arr[i] = (char*)calloc(100, sizeof(char));
    }
    char *code_str = (char*)calloc(100, sizeof(char));
    make_huffman_code(head.pnode, -1, code_str, code_arr, tmp_len);

    write_index_file(data_tmp);

    for(int i = 0; i < MAX_NODE; i++) {
        //printf("%d : %s\n",i,code_arr[i]);
    }

    char *str = (char*)calloc(100, sizeof(char));
    for (int idx = 0; idx < data_tmp.total_freq; idx++) {
        strcpy(str, code_arr[data_tmp.file_data[idx]]);
//        //printf("%s", str);
        len += (int)strlen(str);
        for (int i = 0; str[i]; i++) {
            buffer = buffer << 1;
            buffer = buffer | (unsigned char) (str[i] - '0');
            bit_num += 1;

            if (bit_num == 8) {
                /*
                for(int j=0;j<8;j++) {
                    //printf("%d", (buffer & (1 << j)) >> j);
                }
                //printf("\n");
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
            //printf("%d", (buffer & (1 << j)) >> j);
        }
        //printf("\n");
        */
    }

    //printf("\n");

    for (int i = 0; i < MAX_NODE; i++) {
        free(code_arr[i]);
    }
    free(code_arr);
    free(code_str);

    return len;
}

static void read_fafile(FILE *fp, FILE *w_fp, int tmp_len) {

    const int MAX_NODE = power5(tmp_len);
    char ch_tmp;
    //data info

    Data_info data_tmp; //file data we read
    data_tmp.total_freq = 0;
    data_tmp.data_len = 0;
    data_tmp.data_num = 0;
    data_tmp.freq_arr = (int*)calloc(MAX_NODE + 10, sizeof(int));
    data_tmp.file_data = (int*)calloc((MAX_LENGTH + 10) / tmp_len, sizeof(int));

    int pre_data_num = 0;
    char *odd_data_tmp = (char*)calloc(tmp_len + 2, sizeof(char)); //saving the first char
    odd_data_tmp[0] = '0';
    for (int i = 1; i <= tmp_len + 1; i++) {
        odd_data_tmp[i] = '\0';
    }

    while ((ch_tmp = getc(fp)) != EOF) { //all file reading
        if (ch_tmp == ';' || ch_tmp == '>') { //start of a new data

            if(data_tmp.data_len != 0) { //if we have some read data
                Element head = mk_huffman_tree(data_tmp, odd_data_tmp, tmp_len); //making huffman tree
                //printf("%d end mk huffmantree\n",head.freq);
                encoding(head, data_tmp ,w_fp, tmp_len); //encoding the file
                //printf("%d end encoding\n",data_tmp.data_len);
                destroy_heap(head.pnode);
                //printf("end destroy heap\n");
                free(data_tmp.freq_arr);
                free(data_tmp.file_data);
                free(data_tmp.data_name);
                //printf("end free\n");
                write_index_file(data_tmp);

                //initialization
                data_tmp.data_num = pre_data_num + 1;
                data_tmp.freq_arr = (int*)calloc(MAX_NODE + 10, sizeof(int));
                data_tmp.file_data = (int*)calloc((MAX_LENGTH + 10) / tmp_len, sizeof(int));
                data_tmp.total_freq = 0;
                data_tmp.data_len = 0;
                pre_data_num += 1;
                //printf("end init\n\n\n");
            }
            data_tmp.data_name = read_line(fp, ch_tmp);
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
                if (data_tmp.data_len % tmp_len == 0) { //if it is the first char
                    odd_data_tmp[tmp_len] = ch_tmp;
                    int index = str_to_num(odd_data_tmp + 1, tmp_len);
                    for (int i = 1; i <= tmp_len + 1; i++) {
                        odd_data_tmp[i] = '\0';
                    }
                    data_tmp.freq_arr[index] += 1;
                    data_tmp.file_data[data_tmp.total_freq] = index; //recoding the file info
                    //printf("file_data, idx:%d, %d\n",data_tmp.total_freq,index);
                    data_tmp.total_freq += 1;
                } else { //if it is the second char
                    odd_data_tmp[data_tmp.data_len % tmp_len] = ch_tmp; // record from index 1
                }
            }
        }

        if (data_tmp.data_len == MAX_LENGTH) { //if we read 'MAX_LENGTH' number of char
            Element head = mk_huffman_tree(data_tmp, odd_data_tmp, tmp_len); //make huffman tree
            encoding(head, data_tmp, w_fp, tmp_len); //encoding the file
            //printf("%d end encoding",data_tmp.data_len);
            destroy_heap(head.pnode);
            free(data_tmp.freq_arr);
            free(data_tmp.file_data);
            write_index_file(data_tmp);

            //initialization
            data_tmp.freq_arr = (int*)calloc(MAX_NODE + 10, sizeof(int));
            data_tmp.file_data = (int*)calloc((MAX_LENGTH + 10) / tmp_len, sizeof(int));
            data_tmp.total_freq = 0;
            data_tmp.data_len = 0;
        }
    }

    if ( data_tmp.data_len != 0 ) {
        Element head = mk_huffman_tree(data_tmp, odd_data_tmp, tmp_len); //making huffman tree
        //printf("%d end mk huffmantree\n",head.freq);
        encoding(head, data_tmp, w_fp, tmp_len); //encoding the file
        //printf("%d end encoding\n",data_tmp.data_len);
        destroy_heap(head.pnode);
        //printf("end destroy heap\n");
        free(data_tmp.freq_arr);
        free(data_tmp.file_data);
        free(data_tmp.data_name);
        //printf("end free\n");
        write_index_file(data_tmp);
    }

    free(odd_data_tmp);
}

static void write_index_file(Data_info data_tmp) {

}