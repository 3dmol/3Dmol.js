
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

            <script src="../build/3Dmol.js"></script>
            <style>
             #left{
                float: left;
                width:50%;
            }
            #right{
                float:right;
                width:50%;
            }
            #right ul,#left ul{
                list-style:none;
            }
            .header{
                width:50%;
                position:relative;
                margin: 0 auto;
            }

              #progress {
        position: fixed;
        z-index: 2147483647;
        top: 0;
        left: -6px;
        width: 0%;
        height: 2px;
        background: #b91f1f;
        -moz-border-radius: 1px;
        -webkit-border-radius: 1px;
        border-radius: 1px;
        -moz-transition: width 500ms ease-out,opacity 400ms linear;
        -ms-transition: width 500ms ease-out,opacity 400ms linear;
        -o-transition: width 500ms ease-out,opacity 400ms linear;
        -webkit-transition: width 500ms ease-out,opacity 400ms linear;
        transition: width 500ms ease-out,opacity 400ms linear
    }
    #progress.done {
        opacity: 0
    }
    #progress dd,#progress dt {
        position: absolute;
        top: 0;
        height: 2px;
        -moz-box-shadow: #b91f1f 1px 0 6px 1px;
        -ms-box-shadow: #b91f1f 1px 0 6px 1px;
        -webkit-box-shadow: #b91f1f 1px 0 6px 1px;
        box-shadow: #b91f1f 1px 0 6px 1px;
        -moz-border-radius: 100%;
        -webkit-border-radius: 100%;
        border-radius: 100%
    }
    #progress dd {
        opacity: 1;
        width: 20px;
        right: 0;
        clip: rect(-6px,22px,14px,10px)
    }
    #progress dt {
        opacity: 1;
        width: 180px;
        right: -80px;
        clip: rect(-6px,90px,14px,-6px)
    }
    @-moz-keyframes pulse {
        30% {
            opacity: 1
        }
        60% {
            opacity: 0
        }
        100% {
            opacity: 1
        }
    }
    @-ms-keyframes pulse {
        30% {
            opacity: .6
        }
        60% {
            opacity: 0
        }
        100% {
            opacity: .6
        }
    }
    @-o-keyframes pulse {
        30% {
            opacity: 1
        }
        60% {
            opacity: 0
        }
        100% {
            opacity: 1
        }
    }
    @-webkit-keyframes pulse {
        30% {
            opacity: .6
        }
        60% {
            opacity: 0
        }
        100% {
            opacity: .6
        }
    }
    @keyframes pulse {
        30% {
            opacity: 1
        }
        60% {
            opacity: 0
        }
        100% {
            opacity: 1
        }
    }
    #progress.waiting dd,#progress.waiting dt {
        -moz-animation: pulse 2s ease-out 0s infinite;
        -ms-animation: pulse 2s ease-out 0s infinite;
        -o-animation: pulse 2s ease-out 0s infinite;
        -webkit-animation: pulse 2s ease-out 0s infinite;
        animation: pulse 2s ease-out 0s infinite
    }
        </style>
					</head>
					<body>

          <div id="progress" class="waiting">
    <dt></dt>
    <dd></dd>
    </div>
            <div id="left">
                <ul></ul>
            </div>
            <div id="right">
                <ul></ul>
            </div>
            <div id="gldiv" style="width: 400px; height: 400px; position: relative;"></div>
					"""
end="""<script src="test.js"></script>
	</body>
		 </html>"""
class Example():
    def __init__(self,name,text):
		self.name=name
		self.text=self.parse(text)
    
    def parse(self,text):
        atdata=find_all(text,"@data")
        closecomment=find_all(text,"*/")
        atdiv=find_all(text,"@div")
        if(len(atdata) == 0 and len(atdiv) ==0 ):
            return text

        for data in atdata:
            ending=len(text)
            for com in closecomment:
                if com > data and com<ending:
                    ending=com
            for at in atdiv+atdata:
                if at>data and at<ending:
                    ending=at
            string=text[data+6+len(self.name):ending]
            string="var objectHTML=$(`"+"<textarea style=\\\"display: none;\\\" id=\\\""+self.name+"\\\">"+string+"</textarea>`);document.appendChild(objectHTML);"
            text=text+string

        for data in atdiv:
            ending=len(text)
            for com in closecomment:
                if com > data and com<ending:
                    ending=com
            for at in atdiv+atdata:
                if at>data and at<ending:
                    ending=at
            string=text[data+4:ending]
            string="var objectHTML=$(\\\"`"+string+"`\\\");document.appendChild(objectHTML);"
            text=text+string
        
            print text
        return text



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
	f.write(beggining)
	f.write("<script>system={\n")
	for file in test.files:
		if(type(file.examples)!=type(None) and len(file.examples)>0):
			for example in file.examples:
				f.write(example.name+": function (callback){try{var viewer=$3Dmol.createViewer($(\"#gldiv\"));\n"+example.text+"}catch(err){}},\n")
		elif(type(file.examples)==type(None)):
			f.write(file.contents.name+": function(callback){try{var viewer=$3Dmol.createViewer($(\"#gldiv\"));"+file.contents.text+"}catch(err){}},\n")
	f.write("}</script>")
	f.write(end)



