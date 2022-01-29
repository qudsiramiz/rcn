import os
#folder = os.system("pwd")
folder = "/home/cephadrius/Desktop/git/rcn/figures"
print(folder)
pathiter = (os.path.join(root, filename)
    for root, _, filenames in os.walk(folder)
    for filename in filenames
)

for path in pathiter:
    # the '---' in the example below will be replaced by the '-' in the filenames in the directory
    newname = path.replace('---', '-')
    if newname != path:
        os.rename(path, newname)
        print(path, newname)


