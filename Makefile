# Compiler and flags
CC = cc
CFLAGS = -Werror -Wall -Wextra

# Executable name
NAME = triras

# Source and object files
SRCFILES = tri_raster.c common_func.c
OBJDIR = ./src
OBJFILES = $(SRCFILES:%.c=$(OBJDIR)/%.o)

# Default target
all: $(NAME)

# Link object files to create the final executable
$(NAME): $(OBJFILES)
	$(CC) $(OBJFILES) -o $(NAME)

# Compile source files into object files
$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up generated files
clean:
	rm -f $(OBJFILES)

fclean: clean
	rm -f $(NAME)

# Rebuild everything from scratch
re: fclean all

