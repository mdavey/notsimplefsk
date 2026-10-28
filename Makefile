CC      ?= gcc
CFLAGS  ?= -O2 -Wall -Wextra -std=gnu11 -Isrc
LDLIBS  += -lliquid -lm

BINDIR  := bin
SRC     := src

BINS    := $(BINDIR)/fsk_tx $(BINDIR)/fsk_rx $(BINDIR)/fsk_rx_stream \
           $(BINDIR)/modem_smoke $(BINDIR)/channel

all: $(BINS)

$(BINDIR):
	mkdir -p $(BINDIR)

$(BINDIR)/fsk_tx: $(SRC)/fsk_tx.c $(SRC)/common.c $(SRC)/common.h | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $(SRC)/fsk_tx.c $(SRC)/common.c $(LDLIBS)

$(BINDIR)/fsk_rx: $(SRC)/fsk_rx.c $(SRC)/rx_core.c $(SRC)/rx_core.h \
                   $(SRC)/common.c $(SRC)/common.h | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $(SRC)/fsk_rx.c $(SRC)/rx_core.c $(SRC)/common.c $(LDLIBS)

$(BINDIR)/fsk_rx_stream: $(SRC)/fsk_rx_stream.c $(SRC)/rx_core.c $(SRC)/rx_core.h \
                         $(SRC)/common.c $(SRC)/common.h | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $(SRC)/fsk_rx_stream.c $(SRC)/rx_core.c $(SRC)/common.c $(LDLIBS)

$(BINDIR)/modem_smoke: tests/modem_smoke.c $(SRC)/common.c $(SRC)/common.h | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ tests/modem_smoke.c $(SRC)/common.c $(LDLIBS)

$(BINDIR)/channel: tests/channel.c | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ tests/channel.c $(LDLIBS)

test: all
	sh tests/roundtrip.sh
	sh tests/stream.sh

# manual aplay/arecord loopback (needs a working audio path, not in test)
loopback: all
	sh tests/stream_loopback.sh

# manual loudspeaker -> room -> microphone test (see script header for PHY)
acoustic: all
	sh tests/stream_acoustic.sh

clean:
	rm -rf $(BINDIR)

.PHONY: all test loopback acoustic clean
