CC      ?= gcc
CFLAGS  ?= -O2 -Wall -Wextra -std=gnu11 -Isrc
LDLIBS  += -lliquid -lm

BINDIR  := bin
SRC     := src

BINS    := $(BINDIR)/fsk_tx $(BINDIR)/fsk_rx $(BINDIR)/fsk_rx_stream \
           $(BINDIR)/modem_smoke $(BINDIR)/channel $(BINDIR)/rs_selftest

all: $(BINS)

$(BINDIR):
	mkdir -p $(BINDIR)

$(BINDIR)/fsk_tx: $(SRC)/fsk_tx.c $(SRC)/common.c $(SRC)/common.h \
                   $(SRC)/rs.c $(SRC)/rs.h | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $(SRC)/fsk_tx.c $(SRC)/common.c $(SRC)/rs.c $(LDLIBS)

$(BINDIR)/fsk_rx: $(SRC)/fsk_rx.c $(SRC)/rx_core.c $(SRC)/rx_core.h \
                   $(SRC)/common.c $(SRC)/common.h $(SRC)/rs.c $(SRC)/rs.h | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $(SRC)/fsk_rx.c $(SRC)/rx_core.c $(SRC)/common.c $(SRC)/rs.c $(LDLIBS)

$(BINDIR)/fsk_rx_stream: $(SRC)/fsk_rx_stream.c $(SRC)/rx_core.c $(SRC)/rx_core.h \
                         $(SRC)/common.c $(SRC)/common.h $(SRC)/rs.c $(SRC)/rs.h | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $(SRC)/fsk_rx_stream.c $(SRC)/rx_core.c $(SRC)/common.c $(SRC)/rs.c $(LDLIBS)

$(BINDIR)/modem_smoke: tests/modem_smoke.c $(SRC)/common.c $(SRC)/common.h \
                        $(SRC)/rs.c $(SRC)/rs.h | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ tests/modem_smoke.c $(SRC)/common.c $(SRC)/rs.c $(LDLIBS)

$(BINDIR)/rs_selftest: tests/rs_selftest.c $(SRC)/rs.c $(SRC)/rs.h | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ tests/rs_selftest.c $(SRC)/rs.c

$(BINDIR)/channel: tests/channel.c | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ tests/channel.c $(LDLIBS)

test: all
	$(BINDIR)/rs_selftest
	sh tests/roundtrip.sh
	sh tests/stream.sh
	python3 tests/channel_selftest.py
	sh tests/rx_speed.sh

# realistic impairment/PHY matrix; opt-in because it is much slower than test
stress: all
	python3 tests/stress.py --quick

stress-full: all
	python3 tests/stress.py --full

# full offline RX timing checks, including the 8 KB pathological case
speed-full: all
	RX_SPEED_FULL=1 sh tests/rx_speed.sh

# manual aplay/arecord loopback (needs a working audio path, not in test)
loopback: all
	sh tests/stream_loopback.sh

# manual loudspeaker -> room -> microphone test (see script header for PHY)
acoustic: all
	sh tests/stream_acoustic.sh

clean:
	rm -rf $(BINDIR)

.PHONY: all test stress stress-full speed-full loopback acoustic clean
