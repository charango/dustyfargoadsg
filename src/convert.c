#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>


int main (argc, argv)
int argc;
char *argv[];
{
  struct stat *buf;
  FILE *handle;
  char *array;
  int size, i, j, l;
  char temp;
  if (argc != 2) {
    fprintf (stderr, "Usage : %s filename\n", argv[0]);
    exit (1);
  }
  handle = fopen (argv[1], "r");
  if (handle == NULL) {
    fprintf (stderr, "Can't open %s\n", argv[1]);
    exit (1);
  }
  buf = (struct stat *)malloc(sizeof(struct stat));
  stat (argv[1], buf);
  size = (int)(buf->st_size);
  array = (char *)malloc(size);
  if (array == NULL) {
    fprintf (stderr, "Not enough memory. Sorry.\n", argv[1]);
    exit (1);
  }
  fread (array, sizeof(char), size, handle);
  for (i = 0; i < size/8; i++) {
    for (j = 0; j < 4; j++) {
      l = i*8+j;
      temp = array[l];
      array[l] = array[l+7-2*j];
      array[l+7-2*j] = temp;
    }
  }
  fclose (handle);
  handle = fopen (argv[1], "w");
  if (handle == NULL) {
    fprintf (stderr, "Can't open %s\n", argv[1]);
    exit (1);
  }
  fwrite (array, sizeof(char), size, handle);
  fclose (handle);
  free (array);
  free (buf);
  return 0;
}
  
