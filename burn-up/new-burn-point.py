import datetime

def proc():
    f = open("../rel-todo-gtk3")
    lines = f.readlines()
    add_to_pre_sum = True
    n_above = 0
    n_below = 0
    for line in lines:
        if "------" in line:
            add_to_pre_sum = False
        arr = line.split()
        try:
            if arr[0] == "o":
                if add_to_pre_sum:
                    n_above += 1
                else:
                    n_below += 1
        except IndexError as e:
            pass

    return (n_above, n_below, n_above + n_below)


(n_above, n_below, n_total) = proc()


start = datetime.datetime(2020, 3, 19)
now   = datetime.datetime.now()

delta_s = int(now.strftime("%s")) - int(start.strftime("%s"))
delta_d = delta_s/(60*60*24)
print(delta_d, n_above, n_below, n_total)
