import IPython.display
import time, json

#surface type constants
VDW =1
MS=2
SAS=3
SES=4

class view(object):
    '''A class for constructing embedded 3Dmol.js views in ipython notebooks.
       The results are completely static which means there is no need for there
       to be an active kernel but also that there is no communication between
       the javascript viewer and ipython.
       
       The API for the created object is exactly that for $3Dmol.GLViewer, with
       the exception that the functions all return None.
       http://3dmol.csb.pitt.edu/doc/$3Dmol.GLViewer.html
    '''
    def __init__(self,width=640,height=480,query='',options=dict(),js='http://3dmol.csb.pitt.edu/build/3Dmol.js'):
        '''Create a 3Dmol.js view.
            width -- width in pixels of container
            height -- height in pixels of container
            query -- optional argument to provide to $3Dmol.download
            options -- optional options to provide to $3Dmol.download
            js -- url for 3Dmol.js'''
        divid = "3dmolviewer_UNIQUEID" 
                
        self.startjs = '<div id="%s"  style="position: relative; width: %dpx; height: %dpx">\n' % (divid,width,height)
        self.startjs += '<script>\n'
        self.endjs = '</script>';
        
        #load 3dmol, but only once
        self.startjs += "if(typeof $3Dmolpromise === 'undefined') $3Dmolpromise = $.when($.getScript('%s'))\n" % js
        
        self.startjs += "$3Dmolpromise.done(function() {\n";
        self.endjs = "});\n" + self.endjs
        
        self.startjs += 'var viewer = $3Dmol.createViewer($("#%s"),{backgroundColor:"white"});\n' % divid
        if query:
            self.startjs += '$3Dmol.download("%s", viewer, %s, function() {\n' % (query,json.dumps(options))
            self.endjs = "})\n" + self.endjs        
        self.endjs = "viewer.render();\n" + self.endjs;

    def show(self):
        return IPython.display.HTML(self._repr_html_())
    
    def _repr_html_(self):
        html = (self.startjs+self.endjs).replace('UNIQUEID',str(time.time()).replace('.',''))
        #print html
        return html
    
    def __getattr__(self,name):
        '''auto-instantiate javascript calls based on whatever the user provided'''
        if name.startswith('_'): #object to ipython canary functions
            raise AttributeError("%r object has no attribute %r" %
                         (self.__class__, attr))
        def makejs(*args):            
            cmd = '\tviewer.%s(' % name;
            for arg in args:
                cmd += '%s,' % json.dumps(arg)
            cmd = cmd.rstrip(',')
            cmd += ');\n';
            self.startjs += cmd
            return self
            
        return makejs
