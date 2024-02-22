a = 10 ; b = 3
if ( a == 10 ) and (b == 3) :
    print ('test ok')
s = '''
a = 1
for k in range(10) :
    a = a + 1
print(a)
'''
code = compile(s, '<string>','exec')
exec(code)

name =  input( 'name?' ) 
print (name)
