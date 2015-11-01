str = ''
for i in range(50):
    str += 'a'
for i in range(50):
    str += 'b'
print(str)
step = 20
st = 0
while len(str[st:st+step]) == step:
    for i in range(st):
        print(' ', end = '')
    print(str[st:st+step]) 
    st += 1
    
    
