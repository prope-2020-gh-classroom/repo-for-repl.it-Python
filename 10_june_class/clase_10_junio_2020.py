
string = "xxx Un día soleado como hoy :)"

string2 = "así es!"

print("Seleccionamos string:")
print(string)
print("De la posición 1 a la posición 5:")
print(string[0:5])


string3 = "".join([string, "\n", string2])

print("-"*20)
print("Concatenamos string:")
print(string)
print("con string2:")
print(string2)

print(string3)

print("-"*20)

print("Usamos tres dobles comillas y saltos de línea:")

string4 = """Un día soleado como hoy :)
así es!
"""

print(string4)

print("-"*20)

print("Utilizando splitlines:")
l = string4.splitlines()
print(l)

print("-"*20)

print("Utilizando split:")
print(string4.split())

print("-"*20)

print("Usamos un ciclo for:")
k = 0

for linea in string4.splitlines():
  print(linea)
  print(k)
  k+=1 #k = k + 1

print("-"*20)

print("lista ejemplo:")
l2 = [1,2,3]

for i in l2:
  print(i)



#string5="1 2 3 4 5 6 7 8"
#print(string5.split())







 

