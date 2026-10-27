CC      ?= gcc
CFLAGS  ?= -O2 -Wall -Wextra -std=gnu11 -Isrc
LDLIBS  += -lliquid -lm

BINDIR  := bin
SRC     := src

BINS    := $(BINDIR)/fsk_tx $(BINDIR)/fsk_rx $(BINDIR)/modem_smoke $(BINDIR)/channel

all: $(BINS)

$(BINDIR):
	mkdir -p $(BINDIR)

$(BINDIR)/fsk_tx: $(SRC)/fsk_tx.c $(SRC)/common.c $(SRC)/common.h | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $(SRC)/fsk_tx.c $(SRC)/common.c $(LDLIBS)

$(BINDIR)/fsk_rx: $(SRC)/fsk_rx.c $(SRC)/common.c $(SRC)/common.h | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $(SRC)/fsk_rx.c $(SRC)/common.c $(LDLIBS)

$(BINDIR)/modem_smoke: tests/modem_smoke.c $(SRC)/common.c $(SRC)/common.h | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ tests/modem_smoke.c $(SRC)/common.c $(LDLIBS)

$(BINDIR)/channel: tests/channel.c | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ tests/channel.c $(LDLIBS)

test: all
	sh tests/roundtrip.sh

clean:
	rm -rf $(BINDIR)

.PHONY: all test clean
