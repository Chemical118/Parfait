#include <stdio.h>
#include <stdbool.h>
#include <string.h>

typedef struct Data_info{
    int data_num;
    char *data_name;
    int *freq_arr;
    int total_freq;
    int data_len;
}Data_info;

bool isfafile(char *File); //cheak the file is fasta
char* read_line(FILE *fp); //read one line of the file and return it as a string
void read_fafile(FILE *fp); //read the fasta file
int tmp_to_num(char tmp); // ATGCN to 01234
void mk_huffman_tree(Data_info data_tmp, char odd_data_tmp);

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

int tmp_to_num(char tmp) {
    if (tmp == 'A' || tmp == 'a') {
        return 0;
    } else if (tmp == 'T' || tmp == 't') {
        return 1;
    } else if (tmp == 'G' || tmp == 'g') {
        return 2;
    } else if (tmp == 'C' || tmp == 'c') {
        return 3;
    } else {
        return 4;
    }
}

void read_fafile(FILE *fp) {

    const int MAX_LENGTH = 100000;
    char ch_tmp;
    //data info

    Data_info data_tmp;
    int pre_data_num = 0;
    char odd_data_tmp = '\0';


    while ( (ch_tmp = getc(fp)) != EOF) { //all file reading
        if (ch_tmp == ';' || ch_tmp == '>') { //start of a new data
            mk_huffman_tree(data_tmp, odd_data_tmp);

            int _freq_arr[25] = {};
            data_tmp.data_num = pre_data_num + 1;
            data_tmp.data_name = read_line(fp);
            data_tmp.freq_arr = _freq_arr;
            data_tmp.total_freq = 0;
            data_tmp.data_len = 0;
            pre_data_num += 1;
            if (ch_tmp == ';') {
                read_line(fp);
            }
        } else {
            if (ch_tmp == '\n') continue;
            else {
                data_tmp.data_len += 1;
                if (data_tmp.data_len % 2 == 1) {
                    odd_data_tmp = ch_tmp;
                } else {
                    odd_data_tmp = '\0';
                    data_tmp.freq_arr[tmp_to_num(odd_data_tmp)*5+tmp_to_num(ch_tmp)] += 1;
                    data_tmp.total_freq += 1;
                }
            }
        }

        if (data_tmp.data_len == MAX_LENGTH) {
            mk_huffman_tree(data_tmp, odd_data_tmp);

            int _freq_arr[25] = {};
            data_tmp.freq_arr = _freq_arr;
            data_tmp.total_freq = 0;
            data_tmp.data_len = 0;
        }
    }
}