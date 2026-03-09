FLAGS = -isystem /usr/lib/x86_64-linux-gnu/openmpi/include -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format 
FLAGS_UNUSED = -isystem /usr/lib/x86_64-linux-gnu/openmpi/include -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
 
N = 4
M = 1
P = 1
R = 10
S = 1
TEST_FILES = /home/semyon/BOGACHEVTASK/5_sem/matrix_test/*.txt

all: a.out

a.out: main.cpp f.h f.cpp solve.h solve.cpp io_status.h
	mpicxx $(FLAGS) $^ -o a.out

unused: main.cpp f.h f.cpp solve.h solve.cpp io_status.h
	mpicxx $(FLAGS_UNUSED) $^ -o a.out

zip:
	zip Dubkov_SA.zip *.h *.cpp Makefile -x "test.cpp"


test: a.out
	@echo "Running tests..."
	@for file in $(TEST_FILES); do \
		echo "Testing with $$file:"; \
		mpirun -np $(P) ./a.out $(N) $(M) $(R) 0 $$file; \
		echo "---"; \
	done

test1: a.out
	@echo "Running tests..."
	@for n in $$(seq 1 100); do \
		mpirun -np $(P) ./a.out $$n $(M) $(R) $(S); \
		echo "---"; \
	done

.PHONY: test test1 clean