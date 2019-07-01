#!src/venv/bin/python


from __future__ import print_function
import Tkinter as tk
import sys


class GenericGUI(object):

    def __init__(self, title):
        self.__root = tk.Tk()
        self.__root.winfo_toplevel().title(title)
        self.__initMenuBar()
        self.__initStatusBar()
        frame = tk.Frame(self.__root, bd=1, relief="raised")
        frame.pack(fill="both", expand=1)
        self.__initMainFrame(frame)
        self.__initDisplay(frame)

    def __initMenuBar(self):
        self.__menu_bar = tk.Menu(self.__root)
        self.__root.config(menu=self.__menu_bar)

    def __initStatusBar(self):
        self.__active_status = False
        self.__display_status = tk.StringVar()
        self.toggleActiveStatus(active=False)
        tk.Label(
            self.__root,
            textvariable=self.__display_status,
            relief="sunken",
            anchor="w",
        ).pack(
            side="bottom",
            fill="x"
        )

    def __initMainFrame(self, frame):
        self.main_frame = tk.Frame(
            frame,
        )
        self.main_frame.pack(
            side="left",
            fill="both",
            expand=1
        )

    def __initDisplay(self, frame):
        self.__display_frame = tk.Frame(
            frame,
        )
        self.__display_text = tk.Text(self.__display_frame)
        self.__display_text.config(state="disabled")
        self.__display_text.pack(
            side="left",
            fill="both",
            expand=1
        )
        display_scroll_bar = tk.Scrollbar(
            self.__display_frame,
            orient="vertical"
        )
        display_scroll_bar['command'] = self.__display_text.yview
        display_scroll_bar.pack(
            side="right",
            fill="y"
        )
        self.__display_text['yscrollcommand'] = display_scroll_bar.set
        self.toggleActiveDisplay(active=False)

    def setEasterEggs(self):
        try:
            import src.gui.EasterEgg as easter
            easter.start(self)
        except ImportError:
            pass

    def newSubMenu(self, label):
        sub_menu = tk.Menu(self.__menu_bar, tearoff=0, name=label.lower())
        self.__menu_bar.add_cascade(
            label=label,
            menu=sub_menu
        )
        return sub_menu

    def toggleActiveStatus(self, active=None):
        if active is not None:
            self.__active_status = active
        else:
            self.__active_status = not self.__active_status
        if self.__active_status:
            self.__display_status.set("Busy")
        else:
            self.__display_status.set("Idle")

    def toggleActiveDisplay(self, active=None):
        if active is not None:
            self.__active_display = active
        else:
            self.__active_display = not self.__active_display
        if self.__active_display:
            self.__display_frame.pack(
                side="right",
                fill="both",
            )
        else:
            self.__display_frame.pack_forget()

    def actOnFunction(self, function, *args, **kwargs):
        self.toggleActiveStatus(active=True)
        self.__performFunction(function, *args, **kwargs)
        self.toggleActiveStatus(active=False)

    def actOnFunctions(self, *functions):
        self.toggleActiveStatus(active=True)
        for function in functions:
            self.__performFunction(*function)
        self.toggleActiveStatus(active=False)

    def __performFunction(self, function, *args, **kwargs):
        function(*args, **kwargs)

    def writeTo(self, log_file):
        def write(text):
            self.__display_text.config(state="normal")
            self.__display_text.insert("end", text)
            self.__display_text.update()
            self.__display_text.config(state="disabled")
            if write_options["log_file"] is not None:
                with open(write_options["log_file"], "ab") as log_file:
                    log_file.write(text)
        write_options = {
            "log_file": log_file
        }
        if write_options["log_file"] is not None:
            open(write_options["log_file"], 'wb').close()
        return write

    def start(self):
        self.write = self.writeTo(log_file=None)
        sys.stderr = self
        sys.stdout = self
        self.__root.mainloop()

    def quit(self):
        self.__root.quit()
