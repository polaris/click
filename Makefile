CFLAGS := -O3 -Wall -Weffc++ -pedantic -pedantic-errors -Wextra -Wcast-align \
    -Wcast-qual -Wchar-subscripts -Wcomment -Wconversion -Wdisabled-optimization \
	-Wfloat-equal -Wformat -Wformat=2 -Wformat-nonliteral -Wformat-security -Wformat-y2k \
	-Wimport -Winit-self -Winvalid-pch -Wmissing-braces -Wmissing-field-initializers \
	-Wmissing-format-attribute -Wmissing-include-dirs -Wmissing-noreturn -Wparentheses -Wpointer-arith \
	-Wredundant-decls -Wreturn-type -Wsequence-point -Wshadow -Wsign-compare -Wstack-protector \
	-Wswitch -Wswitch-default -Wswitch-enum -Wtrigraphs -Wuninitialized -Wunknown-pragmas \
	-Wunreachable-code -Wunused -Wunused-function -Wunused-label -Wunused-parameter -Wunused-value \
	-Wunused-variable -Wvariadic-macros -Wvolatile-register-var -Wwrite-strings

LDFLAGS := -lboost_system -lboost_program_options -lasound -lm -lstdc++ -pthread -std=c++14

all: click

click: main.cpp
	$(CC) $(CFLAGS) main.cpp $(LDFLAGS) -o $@

.PHONY: clean
clean:
	@rm -rf *.o click
