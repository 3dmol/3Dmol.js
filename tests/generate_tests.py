
#this program will parse through all of the tests as well as all of the examples and create an html page from that
#the page will have referance images on the right corresponding to the canvas snapshots on the left
#when you clock on the images they will create a webgl instance
#the user should be able to select exactly what tests they wish to run (this can get pretty sophisticated)
#at some point there will be image comparison
import os

def find_all(text,sub):
	examples=[]
	index=0
	while index < len(text):
		index = text.find(sub, index)
		if index == -1:
			break
		examples.append(index)
		index += 2
	return examples
def script_string(typedef,text):
	string="<script>function "+typedef+"{"+text+"}</script>\n"
	return string

beggining="""<!DOCTYPE html>
				<html>
					<head>
						<title></title>
					</head>
					<body>

					"""
end="""</body>
		 </html>"""
class Example():
	def __init__(self,name,text):
		self.name=name
		self.text=text


class File():
	def __init__(self,filename,filetype,contents):
		self.filename=filename
		self.filetype=filetype
		self.contents=Example(filename[4:-3],contents)
		self.examples=None
		if(self.filetype=="generated"):
			self.examples=self.getExamples()

	def getExamples(self):
		nearests=[]
		filename=self.filename
		file=open(filename,'r')
		text=file.read()
		example_indices=find_all(text,"@example")
		typedef_decorators=find_all(text,"@typedef")
		function_decorators=find_all(text,"@function")
		examples=[]
		#name of nearest function or typdef
		name=""
		for i in example_indices:
			#get the index of the next '*/' or the index of the next '@' symbol
			comments=find_all(text,"*/")

			ats=find_all(text,"@")
			smallest_next=0
			for comment in comments:
				if comment>i:
					smallest_next=comment
					break
			for at in ats:
				if(at>i and at<smallest_next):
					smallest_next=at
					break
			#at this point smallest next should be correct
			nearest_dec=-1
			for dec in function_decorators:
				if dec<i and dec>nearest_dec:
					nearest_dec=dec
			for dec in typedef_decorators:
				if dec<i and dec>nearest_dec:
					nearest_dec=dec
			nearests.append(nearest_dec)
			endline=text[nearest_dec:].find("\n")
			if(nearest_dec!=-1):
				name=text[nearest_dec+9:endline+nearest_dec]
				name=name[name.find("#")+1:]
			else:
				name=""
			print name

			exmp=text[i+8:smallest_next]
			exmp=exmp.replace('*','')
			exmp=Example(name,exmp)
			examples.append(exmp)
			file.close()
		count=0
		for i,example in enumerate(examples):
			if(i==0):
				continue
			if(nearests[i]==nearests[i-1]):
				count+=1
				examples[i].name=name+str(count)
			else:
				count=0

		return examples


class TestSystem():
	def __init__(self):
		self.files=self.declareFiles()

	def declareFiles(self):
		manual_tests_path="new/"
		examples_path="../3Dmol/"
		files=[]
		#these are the files with examples in them
		for filename in os.listdir(examples_path):
			if(filename.endswith(".js")):
				text=open(examples_path+filename,'r')
				files.append(File(examples_path+filename,"generated",text.read()))
			elif(filename =="WebGL"):
				path=examples_path+"WebGL/"
				for file in os.listdir(path):
					#paths for all of the files inside of WebGL
					parsed=File(path+file,"generated",open(path+file,"r").read())
					files.append(parsed)
		#these are the build in tests
		for filename in os.listdir(manual_tests_path):
			file=open(manual_tests_path+filename,"r")
			files.append(File(manual_tests_path+filename,"builtin",file.read()))
		
		return files

test=TestSystem()

f=open("one_page.html","w")
f.write("")
f.close()
with open("one_page.html","a") as f:
	for file in test.files:
		if(type(file.examples)!=type(None) and len(file.examples)>0):
			for example in file.examples:
				f.write("<script> function "+example.name+"(){"+example.text+"}</script>\n")
		elif(type(file.examples)==type(None)):
			f.write("<script> function "+file.contents.name+"(){"+file.contents.text+"}</script>\n")




