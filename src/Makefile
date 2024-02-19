CC_FLAGS = gcc -std=c++17 -Wall -Werror -Wextra -g

OS = $(shell uname)

ifeq ($(OS), Linux)
	LINKER_FLAGS = -lgtest -lgtest_main -lrt -lm -lstdc++ -pthread -fprofile-arcs -ftest-coverage
else 
	LINKER_FLAGS = -lgtest -lgtest_main -lm -lpthread -fprofile-arcs -ftest-coverage -lstdc++
endif

all : test

clean : 
	rm -rf *.o *.a *.g* *.out *.info *.dSYM report

rebuild : clean all

test : s21_matrix_oop.a
	$(CC_FLAGS) tests.cpp s21_matrix_oop.a $(LINKER_FLAGS) -L. -ls21_matrix_oop
	./a.out

s21_matrix_oop.a : s21_matrix_oop.o
	ar rcs libs21_matrix_oop.a *.o
	ranlib libs21_matrix_oop.a 
	cp libs21_matrix_oop.a s21_matrix_oop.a

s21_matrix_oop.o : 
	$(CC_FLAGS) -c s21_matrix_oop.cpp

style:
	clang-format -style=google -i *.cpp *.h

style_check:
	clang-format -style=google -n *.cpp *.h

test_leaks: test
	leaks --atExit -- ./a.out

valgrind: test
	 valgrind --vgdb=no --leak-check=full --show-leak-kinds=all --track-origins=yes -s ./a.out
	 make clean

gcov_report : s21_matrix_oop.a
	$(CC_FLAGS) s21*.cpp tests.cpp s21_matrix_oop.a $(LINKER_FLAGS) -L. --coverage -ls21_matrix_oop 
	./a.out
	lcov -t a.out -o report.info --no-external -c -d .
	genhtml -o report report.info
	ifeq ($(OS), Linux)
		xdg-open ./report/index.html
	else
		open ./report/index.html
	endif
