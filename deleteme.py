# delete this

def printme(*args):
    if args:
        print('argument given')
        print(f'args[0] is {args[0]}')
        print(f'type(args) is {type(args)}')
    else:
        print('hello! no argument given')
    return None

# printme(True)

a = 45 
b = 21

def calculateme():
    print(a + b)
    return None

calculateme()