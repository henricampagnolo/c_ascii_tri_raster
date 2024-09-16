CC = cc -Werror -Wall -Wextra
NAME = triras
SRCFILES = tri_raster.c common_func.c
OBJS = ./src
OBJFILES = $(OBJS)/$(SRCFILES).c=.o

all : $(NAME)

$(NAME) : $(OBJFILES)
	$(CC) $(OBJFILES) -o $(NAME)

$(OBJS)/%.o : %.c
	$(CC) -c -o $@ $< 
