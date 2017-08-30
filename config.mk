# Uncomment to statically link with musl
#CC      = musl-gcc
#LDFLAGS = -lm -static -s
#CFLAGS  = -D_GNU_SOURCE -std=c11 -pedantic -Ofast \

CC      = gcc
LDFLAGS = -lm -s
CFLAGS  = -D_GNU_SOURCE -std=c11 -pedantic -O3

#LDFLAGS = -lm
#CFLAGS = -D_GNU_SOURCE -g -static -std=c11 -pedantic -Wall -Wextra -Wunused-macros \
	-Wno-missing-braces -Wno-missing-field-initializers -Wformat=2 \
	-Wswitch-default -Wswitch-enum -Wcast-align -Wpointer-arith \
	-Wbad-function-cast -Wstrict-overflow=5 -Wstrict-prototypes -Winline \
	-Wundef -Wnested-externs -Wcast-qual -Wshadow -Wunreachable-code \
	-Wlogical-op -Wfloat-equal -Wstrict-aliasing=2 -Wredundant-decls \
	-Wold-style-definition -Wmissing-prototypes -Wuninitialized -Wmaybe-uninitialized \
	-Werror -O0 \

