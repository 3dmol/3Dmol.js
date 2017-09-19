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
       
       Optionally, a viewergrid tuple (rows,columns) can be passed to create
       a grid of viewers in a single canvas object.  Successive commands than need to
       specify which viewer they apply to (with viewer=(r,c)) or will apply to all
       viewers in the grid.
       
       The API for the created object is exactly that for $3Dmol.GLViewer, with
       the exception that the functions all return None.
       http://3dmol.csb.pitt.edu/doc/$3Dmol.GLViewer.html
    '''
    def __init__(self,width=640,height=480,query='',viewergrid=None,linked=True,options=dict(),js='https://3dmol.csb.pitt.edu/build/3Dmol.js'):
        '''Create a 3Dmol.js view.
            width -- width in pixels of container
            height -- height in pixels of container
            query -- optional argument to provide to $3Dmol.download
            options -- optional options to provide to $3Dmol.download
            js -- url for 3Dmol.js'''
        divid = "3dmolviewer_UNIQUEID" 
        self.uniqueid = None
        self.startjs = '<div id="%s"  style="position: relative; width: %dpx; height: %dpx">\n' % (divid,width,height)
        self.startjs += '<script>\n'
        self.endjs = '</script>'
        
        self.updatejs = '' # code added since last show
        #load 3dmol, but only once
        self.startjs += "if(typeof $3Dmolpromise === 'undefined') $3Dmolpromise = $.when($.getScript('%s'))\n" % js
        self.startjs += "var viewer_UNIQUEID = null;\n";
        self.startjs += "$3Dmolpromise.done(function() {\n";
        self.endjs = "});\n" + self.endjs
        self.viewergrid = None
        
        if viewergrid: #create split view
            if len(viewergrid) != 2:
                raise ValueError("Incorrectly formated viewergrid arguments.  Must specify rows x columns",viewergrid)
            self.startjs += "var viewergrid_UNIQUEID = null;\n";
            self.startjs += 'viewergrid_UNIQUEID = $3Dmol.createViewerGrid($("#%s"),{rows: %d, cols: %d, control_all: %s},{backgroundColor:"white"});\n' % (divid, viewergrid[0],viewergrid[1],'true' if linked else 'false')
            self.startjs += "viewer_UNIQUEID = viewergrid_UNIQUEID[0][0];\n" 
            self.viewergrid = viewergrid
        else:
            self.startjs += 'viewer_UNIQUEID = $3Dmol.createViewer($("#%s"),{backgroundColor:"white"});\n' % divid
        if query:
            if viewergrid:
                for r in range(viewergrid[0]):
                    for c in range(viewergrid[1]):
                        self.startjs += '$3Dmol.download("%s", viewergrid_UNIQUEID[%d][%d], %s, function() {\n' % (query,r,c,json.dumps(options))
                        self.endjs = "})\n" + self.endjs
            else:
                self.startjs += '$3Dmol.download("%s", viewer_UNIQUEID, %s, function() {\n' % (query,json.dumps(options))
                self.endjs = "})\n" + self.endjs        
        
        if viewergrid:
            for r in range(viewergrid[0]):
                for c in range(viewergrid[1]):
                    self.endjs = "viewergrid_UNIQUEID[%d][%d].render();\n"%(r,c) + self.endjs;
        else:
            self.endjs = "viewer_UNIQUEID.render();\n" + self.endjs;

    def show(self):
        '''Instantiate a new viewer window. Calling this will orphan any previously instantiated viewer windows.'''
        self.updatejs = ''
        return IPython.display.HTML(self._repr_html_())
        
    def insert(self, containerid):
        '''Instead of inserting into notebook here, insert html
        into existing container'''
        html = self._repr_html_()
        html += '''<script>$("#%s").append($("#3dmolviewer_%s")); </script>'''%(containerid,self.uniqueid)
        return IPython.display.HTML(html)
        
    def _repr_html_(self):
        self.uniqueid = str(time.time()).replace('.','')
        self.updatejs = ''
        html = (self.startjs+self.endjs).replace('UNIQUEID',self.uniqueid)
        return html    
        
    def update(self):
        '''Apply commands to existing viewer (will auto-instantiate if necessary).'''
        script = ''
        if self.uniqueid == None:
            script = self._repr_html_()
        script += '''<script>
            $3Dmolpromise.done(function() { //wrap in promise for non-interactive functionality
                %s
                viewer_%s.render();
            });
            </script>''' % (self.updatejs.replace('UNIQUEID',self.uniqueid),self.uniqueid)
        self.updatejs = ''
        return IPython.display.HTML(script)

    def png(self):
        '''output png image of viewer, which must already be instantiated'''
        if not self.uniqueid:
            raise AssertionError('Must instantiate viewer before generating image.')
        script = '''<img id="img_{0}">
            <script>
            var png = $('canvas','#3dmolviewer_{0}')[0].toDataURL();
            $('#img_{0}').attr('src', png)
            </script>'''.format(self.uniqueid)
        return IPython.display.HTML(script)
            
    
    def __getattr__(self,name):
        '''auto-instantiate javascript calls based on whatever the user provided'''
        if name.startswith('_'): #object to ipython canary functions
            raise AttributeError("%r object has no attribute %r" %  (self.__class__, name))
                                                  
        def makejs(*args,**kwargs):            
            if self.viewergrid:
                if 'viewer' in kwargs:
                    coords = kwargs['viewer']
                    if len(coords) != 2 or coords[0] >= self.viewergrid[0] or coords[1] >= self.viewergrid[1] or coords[0] < 0 or coords[1] < 0:
                        raise ValueError("Incorrectly formated viewer argument.  Must specify row and column",self.viewergrid)
                    cmd = '\tviewergrid_UNIQUEID[%d][%d].%s(' % (coords[0],coords[1],name);
                    for arg in args:
                        cmd += '%s,' % json.dumps(arg)
                    cmd = cmd.rstrip(',')
                    cmd += ');\n';
                else: #apply to every viewer
                    cmd = ''
                    for r in range(self.viewergrid[0]):
                        for c in range(self.viewergrid[1]):
                            cmd += '\tviewergrid_UNIQUEID[%d][%d].%s(' % (r,c,name);
                            for arg in args:
                                cmd += '%s,' % json.dumps(arg)
                            cmd = cmd.rstrip(',')
                            cmd += ');\n';
            else:
                cmd = '\tviewer_UNIQUEID.%s(' % name;
                for arg in args:
                    cmd += '%s,' % json.dumps(arg)
                cmd = cmd.rstrip(',')
                cmd += ');\n';

            self.startjs += cmd
            self.updatejs += cmd
            return self
            
        return makejs
