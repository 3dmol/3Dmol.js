#!/usr/bin/env python

#this program will parse through all of the tests as well as all of the examples and create an html page from that
#the page will have referance images on the right corresponding to the canvas snapshots on the left
#when you click on the images they will create a webgl instance
#the user should be able to select exactly what tests they wish to run (this can get pretty sophisticated)
#at some point there will be image comparison
import os, re, glob, sys
#from IPython.core.magics import script
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
        if self.name.startswith('3'): self.name = '_'+self.name
        self.text=self.parse(text)
    
    def parse(self,text):
        self.script = re.sub(re.compile(r'/\*.*?\*/',re.DOTALL),'',text) #example without comments
        self.divs = re.findall(re.compile(r'/\*\s*@div\s*(.*?)\*/',re.DOTALL),text) #all div code
        self.datas = re.findall(re.compile(r'/\*\s*@data\s*(.*?)\*/',re.DOTALL),text) #all data
        self.prescripts = re.findall(re.compile(r'/\*\s*@script\s*(.*?)\*/',re.DOTALL),text) #script to run before div

        #make sure callback is called
        self.script = self.script.replace('viewer.render()','viewer.render(callback)')
        if 'viewer.render(callback)' not in self.script and len(self.divs) == 0:
            self.script += "\nviewer.render(callback);\n";
        
        if self.script.count("viewer.render(callback)") > 1:
            sys.stderr.write("More than one 'viewer.render' call in %s.  FIX THIS NOW!\n"%self.name)
        #construct all javascript test    
        text = ''

        for data in self.datas:
            text += "var wrapper=$.parseHTML(`"+data+"`);\n$('body').append(wrapper);\n"
            
        for div in self.divs:
            text += "var objectHTML=$.parseHTML(`"+div+"`);\n$(\"body\").append(objectHTML);\n"
        if self.divs:
            # remove the div when done with it
            text += "var finished_with_div_callback = function() { callback(viewer); "
            if self.datas:
                text += "$(wrapper).remove(); "
            text += "$(objectHTML).remove();}\n"
            prescript = ''
            for s in self.prescripts:
                prescript += s.replace('viewer.render()','viewer.render(callback)')
            
            if 'viewer.render(callback)' in prescript: # presumably callback in prescript
                prescript = prescript.replace('viewer.render(callback)','viewer.render(finished_with_div_callback)')
                text += prescript
                text += "$3Dmol.autoload(viewer);\n"                
            elif 'viewer.render(callback)' in self.script: # test code calls callback
                self.script.replace('viewer.render(callback)','viewer.render(finished_with_div_callback)')
                text += "$3Dmol.autoload(viewer);\n"
            else:
                text += "$3Dmol.autoload(viewer,finished_with_div_callback);\n"
        
        #code should happen after data is initialized
        text += self.script
        
        return text
        
class File():
    def __init__(self,filename,filetype,contents):
        self.filename=filename
        self.filetype=filetype
        if(self.filetype=="generated"):
            self.examples=self.getExamples()
        else:
            self.examples = [Example(os.path.basename(filename)[0:-3],contents)]

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
            filename=filename.replace("3Dmol/","").lstrip('.').lstrip('/')
            flname=filename+"_"+name
            
            flname=flname.replace(".","_")
            flname=flname.replace("/","_")
            flname=flname.replace("3Dmol_","")
            flname=flname.replace(" ","")
            flname=flname.replace("3Dmol","")
            flname=flname.replace("3","_3")
            flname=flname.replace('$','') #special character in jquery
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
    def __init__(self,d='.'):
        self.files=self.declareFiles(d)

    def declareFiles(self,d):
        manual_tests_path=d+"/tests/auto/tests/"
        examples_path=d+"/3Dmol/"
        files=[]
        #these are the files with examples in them
        for filename in glob.glob(examples_path+'/*.js')+glob.glob(examples_path+'/WebGL/*.js'):
            text=open(filename,'r')
            files.append(File(filename,"generated",text.read()))
        #these are the build in tests
        for filename in glob.glob(manual_tests_path+'/*.js'):
            file=open(filename,"r")
            files.append(File(filename,"builtin",file.read()))
        
        return files

if __name__ == '__main__':
    test=TestSystem() 
    path="tests/auto/build.js"
    
    f=open(path,"w")
    f.write("")
    
    with open(path,"a") as f:
    
        f.write("var system={\n")
    
        for file in test.files:
            for example in file.examples:
                f.write(example.name+": function(viewer,callback,name='"+example.name+"'){try{\n"+example.text+"\n}catch(err){\nconsole.log(err);\nvar li=document.createElement('li');li.innerHTML=name+' '+err;document.getElementById('summary_scroll').appendChild(li);callback(viewer);}},\n")
        f.write("};")



