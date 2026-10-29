"""CustomTkinter GUI for the notsimplefsk G4FSK modem (bin/fsk_tx, bin/fsk_rx_stream).

Run with:  uv run python src/notsimplefsk.py
"""

import os
import queue
import subprocess
import threading
from pathlib import Path
from tkinter import filedialog

import customtkinter as ctk

BIN_DIR = Path(__file__).resolve().parent.parent / "bin"

RATES = ["8000", "16000", "48000", "96000"]

DEFAULTS = {
    "gain": "-12",
    "rate": "48000",
    "fc": "4000",
    "k": "20",
    "h": "0.25",
    "bt": "0.8",
    "delay": "1",
    "preamble": "72",
}

FIELDS = [
    ("gain", "Gain dBFS (TX only)"),
    ("rate", "Sample rate Hz"),
    ("fc", "Carrier freq Hz"),
    ("k", "Samples/symbol"),
    ("h", "Mod index h"),
    ("bt", "BT"),
    ("delay", "Delay (symbols)"),
    ("preamble", "Preamble (symbols)"),
]


def bin_path(name):
    return str(BIN_DIR / name)


def _num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def derived(rate, k, h):
    """Return (baud, bitrate_bps, occupied_bw_hz); None components on bad input.

    Matches the modem's own math (src/common.c): Rs = rate/k, 2 bits/symbol,
    Carson BW = 2 * (1.5*h*Rs + Rs).
    """
    rate, k, h = _num(rate), _num(k), _num(h)
    baud = None
    if rate is not None and k not in (None, 0.0):
        baud = rate / k
    bitrate = 2.0 * baud if baud is not None else None
    bw = 2.0 * (1.5 * h * baud + baud) if baud is not None and h is not None else None
    return baud, bitrate, bw


def build_args(kind, params):
    """Build argv for fsk_tx/fsk_rx_stream from raw UI strings (verbatim)."""
    args = [bin_path("fsk_tx" if kind == "tx" else "fsk_rx_stream")]
    if kind == "tx":
        args += ["-g", params["gain"]]
    args += ["-r", params["rate"], "-c", params["fc"], "-k", params["k"],
             "-H", params["h"], "-b", params["bt"], "-d", params["delay"],
             "-P", params["preamble"]]
    if kind != "tx":
        args.append("--verbose")
    return args


class App(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("notsimplefsk modem")
        self.geometry("560x640")
        self._queue = queue.Queue()
        self._build_ui()
        self._refresh_derived()
        self.after(100, self._poll)

    def _build_ui(self):
        settings = ctk.CTkFrame(self)
        settings.pack(fill="x", padx=12, pady=(12, 6))

        self.entries = {}
        for row, (key, label) in enumerate(FIELDS):
            ctk.CTkLabel(settings, text=label, anchor="e").grid(
                row=row, column=0, sticky="e", padx=(12, 8), pady=4)
            if key == "rate":
                # CTkOptionMenu calls its command with the new value; swallow it.
                w = ctk.CTkOptionMenu(settings, values=RATES, width=140,
                                      command=lambda _v: self._refresh_derived())
                w.set(DEFAULTS[key])
                self.rate_menu = w
            else:
                w = ctk.CTkEntry(settings, width=140)
                w.insert(0, DEFAULTS[key])
                w.bind("<KeyRelease>", lambda e: self._refresh_derived())
                self.entries[key] = w
            w.grid(row=row, column=1, sticky="w", padx=(0, 12), pady=4)

        self.derived_label = ctk.CTkLabel(self, text="", anchor="w")
        self.derived_label.pack(fill="x", padx=16, pady=(2, 6))

        buttons = ctk.CTkFrame(self, fg_color="transparent")
        buttons.pack(fill="x", padx=12)
        self.tx_btn = ctk.CTkButton(buttons, text="TX: data -> audio",
                                    command=lambda: self._start("tx"))
        self.tx_btn.pack(side="left", expand=True, fill="x", padx=(0, 6))
        self.rx_btn = ctk.CTkButton(buttons, text="RX: audio -> data",
                                    command=lambda: self._start("rx"))
        self.rx_btn.pack(side="left", expand=True, fill="x", padx=(6, 0))
        self.play_btn = ctk.CTkButton(buttons, text="Play",
                                      command=self._play)
        self.play_btn.pack(side="left", expand=True, fill="x", padx=(6, 0))

        self.log = ctk.CTkTextbox(self, height=200, state="disabled",
                                  wrap="word", font=("TkFixedFont", 12))
        self.log.pack(fill="both", expand=True, padx=12, pady=12)

    def _refresh_derived(self):
        baud, bitrate, bw = derived(
            self._params()["rate"], self.entries["k"].get(), self.entries["h"].get())
        f = lambda v, unit="", scale=1.0, prec=0: (
            "—" if v is None else f"{v / scale:,.{prec}f}{unit}")
        self.derived_label.configure(
            text=f"Baud: {f(baud)}  |  Bitrate: {f(bitrate, ' bps')}  |  "
                 f"Occupied BW: {f(bw, ' kHz', 1000.0, 2)}")

    def _params(self):
        params = {key: w.get() for key, w in self.entries.items()}
        params["rate"] = self.rate_menu.get()
        return params

    def _start(self, kind):
        s16le = [("s16le audio", "*.s16le")]
        if kind == "tx":
            in_title, in_types = "Select input data file", [("all files", "*.*")]
            out_title, out_ext, out_types = "Save audio output as", ".s16le", s16le
        else:
            in_title, in_types = "Select input audio file", s16le
            out_title, out_ext = "Save recovered data as", ".bin"
            out_types = [("all files", "*.*")]
        in_path = filedialog.askopenfilename(
            title=in_title, filetypes=in_types, parent=self)
        if not in_path:
            return
        out_path = filedialog.asksaveasfilename(
            title=out_title, defaultextension=out_ext, filetypes=out_types,
            parent=self)
        if not out_path:
            return
        args = build_args(kind, self._params())
        self.tx_btn.configure(state="disabled")
        self.rx_btn.configure(state="disabled")
        threading.Thread(target=self._run, args=(args, in_path, out_path),
                         daemon=True).start()

    def _play(self):
        path = filedialog.askopenfilename(
            title="Select audio file to play",
            filetypes=[("s16le audio", "*.s16le")],
            parent=self)
        if not path:
            return
        args = ["ffplay", "-f", "s16le", "-ar", self._params()["rate"], path]
        try:
            proc = subprocess.Popen(args, stdin=subprocess.DEVNULL)
        except FileNotFoundError:
            self._queue.put(("$ " + " ".join(args)
                             + "\nerror: ffplay not found", False))
            return
        # reap the child so repeated Play clicks do not leave zombies
        threading.Thread(target=proc.wait, daemon=True).start()
        self._queue.put(("$ " + " ".join(args) + "\nstarted ffplay", False))

    def _run(self, args, in_path, out_path):
        lines = ["$ " + " ".join(args),
                 f"stdin:  {in_path}",
                 f"stdout: {out_path}"]
        try:
            with open(in_path, "rb") as f:
                data = f.read()
            proc = subprocess.run(args, input=data, capture_output=True)
            with open(out_path, "wb") as f:
                f.write(proc.stdout)
            lines.append(f"return code: {proc.returncode}")
            err = proc.stderr.decode("utf-8", "replace").strip()
            if err:
                lines.append("stderr:\n" + err)
            lines.append(f"wrote {os.path.getsize(out_path)} bytes")
        except FileNotFoundError:
            lines.append(f"error: {args[0]} not found "
                         f"(run `make` in {BIN_DIR.parent})")
        except OSError as e:
            lines.append(f"error: {e}")
        except Exception as e:                      # keep the UI usable
            lines.append(f"error: unexpected {type(e).__name__}: {e}")
        finally:
            self._queue.put(("\n".join(lines), True))

    def _poll(self):
        try:
            while True:
                msg, done = self._queue.get_nowait()
                self.log.configure(state="normal")
                self.log.insert("end", msg + "\n\n")
                self.log.see("end")
                self.log.configure(state="disabled")
                if done:
                    self.tx_btn.configure(state="normal")
                    self.rx_btn.configure(state="normal")
        except queue.Empty:
            pass
        self.after(100, self._poll)


def main():
    App().mainloop()


if __name__ == "__main__":
    main()
