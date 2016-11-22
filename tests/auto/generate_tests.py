
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
def find_next(text,sub,currenti):
    index=currenti
    i=-1
    while index < len(text):
        index = text.find(sub, index)
        if index == -1:
            break
        i=index
        index += 2
    return i
def script_string(typedef,text):
    string="function "+typedef+"{"+text+"}\n"
    return string

class Example():
    def __init__(self,name,text):
        self.name=name.replace(".","_")
        self.name=self.name.replace("/","_")
        self.name=self.name.replace("_3Dmol_","")
        self.name=self.name.replace("___","")
        self.name=self.name.replace("_tests_","")
        self.text=self.parse(text)
    
    def parse(self,text):
        atdata=find_all(text,"@data")
        closecomment=find_all(text,"*/")
        atdiv=find_all(text,"@div")
        text=text.replace("myviewer","viewer")
        if(find_all(text,"viewer.render()")==[] and find_all(text,"viewer.render(callback)")==[] and find_all(text,"@data")==[] and find_all(text,"@div")==[]):
            text=text+"viewer.render(callback);"
        if(len(atdata) == 0 and len(atdiv) ==0 ):
            text=text.replace("viewer.render()","viewer.render(callback)")
            return text

        for data in atdata:
            ending=len(text)
            for com in closecomment:
                if com > data and com<ending:
                    ending=com
            for at in atdiv+atdata:
                if at>data and at<ending:
                    ending=at

            string=text[data+6:ending]
            text=text[0:data]+text[ending:]
            string="var wrapper=$.parseHTML(`"+string+"`);\n$('#div_"+self.name+"').append(wrapper);\nviewer.autoload();";
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
            text=text[0:data]+text[ending:]
            string="var objectHTML=$.parseHTML(`"+string+"`);$(\"body\").append(objectHTML);\nglobal_viewer=viewer;\nglobal_callback=callback;\nviewer.autoload(viewer);"
            text=text+string
        text=text.replace("viewer.render()","viewer.render(callback)")
        return text
        
class File():
    def __init__(self,filename,filetype,contents):
        self.filename=filename
        self.filetype=filetype
        self.contents=Example(filename[10:-3],contents)
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
            filename=filename.replace("3Dmol/","")
            flname=filename+"_"+name
            flname=flname.replace(".","_")
            flname=flname.replace("/","_")
            flname=flname.replace("3Dmol_","")
            flname=flname.replace(" ","")
            flname=flname.replace("3Dmol","")
            flname=flname.replace("3","_3")
            exmp=Example(flname,exmp)
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
        manual_tests_path="tests/auto/tests/"
        examples_path="3Dmol/"
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
        exceptions=["generate_tests.py","tests.html","test.js","volData","imgs","url.cgi","tests.html"]
        for filename in os.listdir(manual_tests_path):
            if(filename in exceptions):
                continue
            file=open(manual_tests_path+filename,"r")
            files.append(File(manual_tests_path+filename,"builtin",file.read()))
        
        return files

test=TestSystem() 
path="tests/auto/build.js"

f=open(path,"w")
f.write("")

with open(path,"a") as f:

    f.write("""var global_viewer=null;
               var global_callback=null;
                    function div_callback(){
                        global_viewer.render(global_callback);

                    }""")
    f.write("var system={\n")

    for file in test.files:
        if(type(file.examples)!=type(None) and len(file.examples)>0):
            for example in file.examples:
                f.write(example.name+": function(viewer,callback,name='"+example.name+"'){try{\n"+example.text+"\n}catch(err){callback();}},\n")
        elif(type(file.examples)==type(None)):
            f.write(file.contents.name+": function(viewer,callback,name='"+file.contents.name+"'){try{\n"+file.contents.text+"\n}catch(err){callback();}},\n")
    f.write("};")



