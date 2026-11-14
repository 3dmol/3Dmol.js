#!/usr/bin/env python3

#this program will parse through all of the tests as well as all of the examples and create an html page from that
#the page will have referance images on the right corresponding to the canvas snapshots on the left
#when you click on the images they will create a webgl instance
#the user should be able to select exactly what tests they wish to run (this can get pretty sophisticated)
#at some point there will be image comparison
import os, re, glob, sys, json, subprocess
from posixpath import abspath
from pathlib import Path
#from IPython.core.magics import script

def script_string(typedef,text):
    string="function "+typedef+"{"+text+"}\n"
    return string

def extract_files(files):
    """Parse all files in one Node process using the installed TypeScript parser."""
    helper = Path(__file__).with_name("extract_examples.cjs")
    result = subprocess.run(["node", str(helper)], input=json.dumps(files),
                            text=True, capture_output=True)
    if result.returncode:
        raise ValueError("Unable to extract render examples:\n" + result.stderr)
    return json.loads(result.stdout)


class Example():
    def __init__(self, name, text, prepared=None):
        self.name=name.replace(".","_")
        self.name=self.name.replace("/","_")
        self.name=self.name.replace("_3Dmol_","")
        self.name=self.name.replace("___","")
        self.name=self.name.replace("_tests_","")
        if self.name.startswith('3'): self.name = '_'+self.name
        if prepared is None:
            prepared = extract_files([{"filename": name + ".js", "kind": "builtin", "text": text}])[0]["examples"][0]
        self.text=self.parse(prepared)

    def parse(self, prepared):
        self.script = prepared["script"]
        self.divs = prepared["divs"]
        self.datas = prepared["datas"]
        self.prescripts = prepared["prescripts"]

        #construct all javascript test
        script = self.script
        text = ''

        if self.datas:
            text += "var dataWrappers=[];\n"
        for data in self.datas:
            text += "var wrapper=$.parseHTML(`"+data+"`);\ndataWrappers.push(wrapper);\n$('body').append(wrapper);\n"
            
        for div in self.divs:
            text += "var objectHTML=$.parseHTML(`"+div+"`);\n$(\"body\").append(objectHTML);\n"
        if self.datas:
            text += "var original_callback=callback;\ncallback=function(viewer){try{original_callback(viewer);}finally{dataWrappers.forEach(function(wrapper){$(wrapper).remove();});}};\n"
        if self.divs:
            # remove the div when done with it
            text += "var original_div_callback = callback;\ncallback = function(viewer) { try { original_div_callback(viewer); } finally { $(objectHTML).remove(); } };\n"
            prescript = ''
            for s in self.prescripts:
                prescript += s
            
            if 'viewer.render(callback)' in prescript: # presumably callback in prescript
                text += prescript
                text += "$3Dmol.autoload(viewer);\n"                
            elif 'viewer.render(callback)' in self.script: # test code calls callback
                text += "$3Dmol.autoload(viewer);\n"
            else:
                text += "$3Dmol.autoload(viewer,callback);\n"
        
        #code should happen after data is initialized
        text += script
        
        return text
        
class File():
    def __init__(self, filename, filetype, contents, extracted=None):
        self.filename = filename
        self.filetype = filetype
        if extracted is None:
            extracted = extract_files([{"filename": filename, "kind": filetype, "text": contents}])[0]["examples"]
        self.examples = self.getExamples(contents, extracted)

    def getExamples(self, text, extracted=None):
        if extracted is None:
            extracted = extract_files([{"filename": self.filename, "kind": self.filetype, "text": text}])[0]["examples"]
        filename = re.sub(r'.*/src/', "", self.filename)
        examples = []
        for entry in extracted:
            if self.filetype == "generated":
                # Preserve the historical first-example names and existing image
                # baselines. Additional examples use their one-based ordinal.
                name = filename + "_" + entry["name"]
                name = name.replace(".", "_").replace("/", "_")
                name = name.replace("3Dmol_", "").replace(" ", "")
                name = name.replace("3Dmol", "").replace("3", "_3").replace('$', '')
                if entry["ordinal"] > 1:
                    name += "_" + str(entry["ordinal"])
            else:
                name = os.path.basename(self.filename)[:-3]
            examples.append(Example(name, entry["code"], entry))
        return examples


class TestSystem():
    def __init__(self, d):
        self.files = self.declareFiles(d)

    def declareFiles(self, d):
        root = Path(d).resolve()
        sources = sorted((root / "src").rglob("*.ts")) + sorted((root / "src").rglob("*.js"))
        manual = sorted((root / "tests/auto/tests").glob("*.js"))
        inputs = []
        for filename in sources + manual:
            text = filename.read_text(encoding="utf-8")
            kind = "builtin" if filename in manual else "generated"
            if kind == "builtin" and '/* skip */' in text:
                continue
            inputs.append({"filename": str(filename), "kind": kind, "text": text})
        extracted = extract_files(inputs)
        files = [File(source["filename"], source["kind"], source["text"], result["examples"])
                 for source, result in zip(inputs, extracted)]
        owners = {}
        for file in files:
            for example in file.examples:
                if example.name in owners:
                    raise ValueError("Duplicate render example name %r in %s and %s" %
                                     (example.name, owners[example.name], file.filename))
                owners[example.name] = file.filename
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
                f.write(example.name+": function(viewer,callback,name='"+example.name+"'){try{\n"+example.text+"\n}catch(err){\nconsole.log(err);\nvar summary=document.getElementById('summary_list');if(summary){var li=document.createElement('li');li.textContent=name+' '+err;summary.appendChild(li);}throw err;}},\n")
        f.write("};")

    # The web server needs only Python and this build artifact, not Node or
    # TypeScript on its CGI PATH. Keep one extractor for every test harness.
    manifest = {
        "version": 1,
        "examples": {
            example.name: {field: getattr(example, field)
                           for field in ("script", "divs", "datas", "prescripts")}
            for file in test.files for example in file.examples
        }
    }
    Path(path).with_name("examples.json").write_text(
        json.dumps(manifest, ensure_ascii=False) + "\n", encoding="utf-8")
