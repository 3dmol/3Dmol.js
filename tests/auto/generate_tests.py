#!/usr/bin/env python3

#this program will parse through all of the tests as well as all of the examples and create an html page from that
#the page will have referance images on the right corresponding to the canvas snapshots on the left
#when you click on the images they will create a webgl instance
#the user should be able to select exactly what tests they wish to run (this can get pretty sophisticated)
#at some point there will be image comparison
import os, re, glob, sys
from posixpath import abspath
#from IPython.core.magics import script

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
        if not re.search(r'viewer.render(.*callback)', self.script) and len(self.divs) == 0:
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
            self.examples=self.getExamples(contents)
        else:
            self.examples = [Example(os.path.basename(filename)[0:-3],contents)]

    def getExamples(self, text):
        nearests=[]
        filename=self.filename
        filename=re.sub(r'.*/src/',"",filename)

        #extract  example block from before function/methods
        fex =  re.findall(r'@example(((?!\*/).)*?)\*/\s+(export\s+function|public|\s*)\s+(\S+)\s*\(',text,re.DOTALL)
         #extract just code and name
        fex = [(m[0],m[-1]) for m in fex]

        #and from typedefs
        tex = re.findall(r'@typedef\s+(\S+)\s+((?!\*/).)*@example(.*?)\*/',text,re.DOTALL)
        tex = [(m[-1],m[0]) for m in tex]
        
        #classes
        cex =  re.findall(r'@example(((?!\*/).)*?)\*/\s+export (class|const|enum|interface)\s+(\S+)',text,re.DOTALL)
        cex = [(m[0],m[-1]) for m in cex]


        def makename(name):
            #massage file/func name
            flname=filename+"_"+name
            
            flname=flname.replace(".","_")
            flname=flname.replace("/","_")
            flname=flname.replace("3Dmol_","")
            flname=flname.replace(" ","")
            flname=flname.replace("3Dmol","")
            flname=flname.replace("3","_3")
            flname=flname.replace('$','') #special character in jquery

            return flname

        examples = []
        for code,name in fex+tex+cex:
            code = code.replace('*','')
            name = makename(name)
            #look for @directives and stop the code there
            m = re.search('(^.*?)@(\S+).*',code,re.DOTALL)
            if m:
                code = m.group(1)
                if m.group(2) == 'example':
                    print(f"Need to handle multiple examples in {name}")
            examples.append(Example(name,code))


        return examples


class TestSystem():
    def __init__(self,d):
        self.files=self.declareFiles(d)

    def declareFiles(self,d):
        manual_tests_path=abspath(d+"/tests/auto/tests/")
        examples_path=abspath(d+"/src")
        files=[]
        # typescript files
        for filename in glob.glob(examples_path+'/*.ts'):
            with open(filename, 'r', encoding='utf-8') as file:
                contents = file.read()
                files.append(File(filename,"generated",contents))

        #these are the files with examples in them
        for filename in glob.glob(examples_path+'/*.js'):
            with open(filename,'r', encoding="utf-8") as text:
                files.append(File(filename,"generated",text.read()))
        #these are the built in tests
        for filename in glob.glob(manual_tests_path+'/*.js'):
            with open(filename,"r", encoding="utf-8") as file:
                files.append(File(filename,"builtin",file.read()))
        
        #print("Found %d files"%len(files))
        
        return files

if __name__ == '__main__':
    
    rootdir = os.getcwd()
    if rootdir.endswith('tests/auto'): # in test dir
        path = abspath(rootdir + '/build.js')
        rootdir = abspath("../..")
    else: # assume run from 3Dmol        
        path=abspath("./tests/auto/build.js")

    test=TestSystem(rootdir) 
    f=open(path,"w")
    f.write("")
    
    with open(path,"a") as f:
    
        f.write("var system={\n")
    
        for file in test.files:
            for example in file.examples:
                f.write(example.name+": function(viewer,callback,name='"+example.name+"'){try{\n"+example.text+"\n}catch(err){\nconsole.log(err);\nvar li=document.createElement('li');li.innerHTML=name+' '+err;document.getElementById('summary_scroll').appendChild(li);callback(viewer);}},\n")
        f.write("};")



