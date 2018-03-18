def get_extension1(filename):
    return(filename.split(".")[-1])

def get_extension2(filename):
    import os.path
    return(os.path.splitext(filename)[1])

def get_extension3(filename):
    return filename[filename.rfind('.'):][1:]

print(get_extension1("hello.txt.gz"))
print(get_extension1("hello"))