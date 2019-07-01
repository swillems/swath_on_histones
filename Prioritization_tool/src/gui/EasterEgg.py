#!src/venv/bin/python


import tkMessageBox


def start(gui):
    def getPaperMenu():
        if "paper" not in gui._GenericGUI__menu_bar.children:
            return gui.newSubMenu(label="Paper")
        return gui._GenericGUI__menu_bar.children["paper"]

    def natureEgg(*args):
        if "natureEgg" in eggs:
            return
        eggs.add("natureEgg")
        paper_menu = getPaperMenu()
        paper_menu.add_command(
            label="Auto-generate Nature Paper",
            command=lambda: tkMessageBox.showinfo(
                "Auto-generate Nature Paper",
                "This will be implemented in DHAENENS 3.0!"
            )
        )

    def scienceEgg(*args):
        if "scienceEgg" in eggs:
            return
        eggs.add("scienceEgg")
        paper_menu = getPaperMenu()
        paper_menu.add_command(
            label="Auto-generate Science Paper",
            command=lambda: tkMessageBox.showinfo(
                "Auto-generate Science Paper",
                "This will be implemented in DHAENENS 3.0!"
            )
        )

    def cellEgg(*args):
        if "cellEgg" in eggs:
            return
        eggs.add("cellEgg")
        paper_menu = getPaperMenu()
        paper_menu.add_command(
            label="Auto-generate Cell Paper",
            command=lambda: tkMessageBox.showinfo(
                "Auto-generate Cell Paper",
                "This will be implemented in DHAENENS 3.0!"
            )
        )

    def jprEgg(*args):
        if "jprEgg" in eggs:
            return
        eggs.add("jprEgg")
        paper_menu = getPaperMenu()
        paper_menu.add_command(
            label="Auto-generate JPR Paper",
            command=lambda: tkMessageBox.showinfo(
                "Auto-generate JPR Paper",
                "This will be implemented in DHAENENS 3.0!"
            )
        )

    def findEgg(*args):
        eggs.add("findEgg")
        tkMessageBox.showinfo(
            "Ctrl-f",
            "Unfortunately no search function is available.\n"
            "However -if you really want to find something-, "
            "there might just be some Easter Eggs to be found...\n"
            "\n"
            "Including this one you have found {}/{} Easter Eggs!".format(
                len(eggs),
                len(eggs_to_find)
            )
        )

    eggs_to_find = [
        ("<Control-f>", findEgg),
        ("<Control-Shift-KeyPress-N>", natureEgg),
        ("<Control-Shift-KeyPress-S>", scienceEgg),
        ("<Control-Shift-KeyPress-C>", cellEgg),
        ("<Control-Shift-KeyPress-J>", jprEgg),
    ]
    eggs = set()
    for egg, egg_function in eggs_to_find:
        gui._GenericGUI__root.bind(egg, egg_function)


# def start(gui):
#     keys = set()
#     def keyPressHandler(event):
#         if event.char not in keys:
#             print(keys)
#         keys.add(event.char)
#         if keys == set("SANDER"):
#             tkMessageBox.showinfo(
#                 "Auto-generate Nature Paper",
#                 "This will be implemented in DHAENENS 3.0!"
#             )
#     def keyReleaseHandler(event):
#         keys.remove(event.char)
#         # keys=set()
#     gui.root.bind_all('<KeyPress>', keyPressHandler)
#     # gui.root.bind_all('<Control-Shift>', keyPressHandler)
#     gui.root.bind_all('<KeyRelease>', keyReleaseHandler)
