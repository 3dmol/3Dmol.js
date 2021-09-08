import time
import json
import sys

try:
    import IPython.display
    _has_IPython = True
except ImportError:
    _has_IPython = False

if sys.version_info >= (3, 8):
    from importlib import metadata
else:
    import importlib_metadata as metadata

__version__ = metadata.version('py3Dmol')

#surface type constants
VDW =1
MS=2
SAS=3
SES=4

def using_ipython(func):

    def inner(*args, **kwargs):
        if not _has_IPython:
            raise ImportError("This function requires an active IPython notebook.")
        return func(*args, **kwargs)

    return inner

def tostr(arg):
    '''convert an argument to javascript string,
    can handle top-level numpy/torch arrays'''
    try:
        return json.dumps(arg)
    except TypeError:
        return json.dumps(arg.tolist()) #numpy arrays and torch tensor
        
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
       http://3dmol.org/doc/$3Dmol.GLViewer.html
    '''
    def __init__(self,width=640,height=480,query='',viewergrid=None,data=None,style=None,linked=True,options=dict(),js='https://3dmol.org/build/3Dmol.js'):
        '''Create a 3Dmol.js view.
            width -- width in pixels of container
            height -- height in pixels of container
            query -- optional argument to provide to $3Dmol.download
            viewergrid -- optional tuple (rows,columns) to define grid
            data -- molecular data to provide to addModel, wit viewer grid can be indexed (r,c)
            style -- style to apply, with viewer grid can be indexed (r,c)
            options -- optional options to provide to $3Dmol.download
            js -- url for 3Dmol.js'''
        divid = "3dmolviewer_UNIQUEID"
        warnid = "3dmolwarning_UNIQUEID"
        self.uniqueid = None
        self.startjs = '''<div id="%s"  style="position: relative; width: %dpx; height: %dpx">
        <p id="%s" style="background-color:#ffcccc;color:black">You appear to be running in JupyterLab (or JavaScript failed to load for some other reason).  You need to install the 3dmol extension: <br>
        <tt>jupyter labextension install jupyterlab_3dmol</tt></p>
        </div>\n''' % (divid,width,height,warnid)
        self.startjs += '<script>\n'
        self.endjs = '</script>'

        self.updatejs = '' # code added since last show
        #load 3dmol, but only once, can't use jquery :-(
        #https://medium.com/@vschroeder/javascript-how-to-execute-code-from-an-asynchronously-loaded-script-although-when-it-is-not-bebcbd6da5ea
        self.startjs += """
var loadScriptAsync = function(uri){
  return new Promise((resolve, reject) => {
    var tag = document.createElement('script');
    tag.src = uri;
    tag.async = true;
    tag.onload = () => {
      resolve();
    };
  var firstScriptTag = document.getElementsByTagName('script')[0];
  firstScriptTag.parentNode.insertBefore(tag, firstScriptTag);
});
};

if(typeof $3Dmolpromise === 'undefined') {
$3Dmolpromise = null;
  $3Dmolpromise = loadScriptAsync('%s');
}

var viewer_UNIQUEID = null;
var warn = document.getElementById("%s");
if(warn) {
    warn.parentNode.removeChild(warn);
}
""" % (js,warnid)

        self.startjs += "$3Dmolpromise.then(function() {\n";
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
                    cmds = ''
                    if data:
                        try:
                            d = data[r][c]
                        except:
                            d = data
                        self.startjs += "viewergrid_UNIQUEID[%d][%d].addModel(%s);\n"%(r,c,json.dumps(d))
                    if style:
                        try:
                            s = style[r][c]
                        except:
                            s = style
                        cmds += "viewergrid_UNIQUEID[%d][%d].setStyle(%s);\n"%(r,c,json.dumps(s))
                    self.startjs += cmds+"viewergrid_UNIQUEID[%d][%d].zoomTo();"%(r,c) 
                    self.endjs = "viewergrid_UNIQUEID[%d][%d].render();\n"%(r,c)+self.endjs;
        else:
            cmds = ''
            if data:
                cmds = "viewer_UNIQUEID.addModel(%s);\n"%json.dumps(data)
            if style:
                cmds += "viewer_UNIQUEID.setStyle(%s);\n"%json.dumps(style)
            self.startjs += cmds + "viewer_UNIQUEID.zoomTo();\n"
            self.endjs = "viewer_UNIQUEID.render();\n" + self.endjs;
            

    @using_ipython
    def show(self):
        '''Instantiate a new viewer window. Calling this will orphan any previously instantiated viewer windows.'''
        self.updatejs = ''
        html = self._make_html()
        return IPython.display.publish_display_data({'application/3dmoljs_load.v0':html, 'text/html': html},metadata={})

    @using_ipython
    def insert(self, containerid):
        '''Instead of inserting into notebook here, insert html
        into existing container'''
        html = self._make_html()
        html += '''<script>$("#%s").append($("#3dmolviewer_%s")); </script>'''%(containerid,self.uniqueid)
        return IPython.display.publish_display_data({'application/3dmoljs_load.v0':html, 'text/html': html},metadata={})

    def _make_html(self):
        self.uniqueid = str(time.time()).replace('.','')
        self.updatejs = ''
        html = (self.startjs+self.endjs).replace('UNIQUEID',self.uniqueid)
        return html

    @using_ipython
    def _repr_html_(self):
        html = self._make_html()
        return IPython.display.publish_display_data({'application/3dmoljs_load.v0':html, 'text/html': html},metadata={})

    @using_ipython
    def update(self):
        '''Apply commands to existing viewer (will auto-instantiate if necessary).'''
        script = ''
        if self.uniqueid == None:
            script = self._make_html()
        script += '''<script>
            $3Dmolpromise.then(function() { //wrap in promise for non-interactive functionality
                %s
                viewer_%s.render();
            });
            </script>''' % (self.updatejs.replace('UNIQUEID',self.uniqueid),self.uniqueid)
        self.updatejs = ''
        return IPython.display.publish_display_data({'application/3dmoljs_load.v0':script, 'text/html': script},metadata={})

    @using_ipython
    def png(self):
        '''output png image of viewer, which must already be instantiated'''
        if not self.uniqueid:
            raise AssertionError('Must instantiate viewer before generating image.')
        script = '''<img id="img_{0}">
            <script>
            var png = viewer_{0}.pngURI()
            $('#img_{0}').attr('src', png)
            </script>'''.format(self.uniqueid)
        return IPython.display.publish_display_data({'application/3dmoljs_load.v0':script, 'text/html': script},metadata={})
        
    class model(object):
      '''Wrapper for referencing a model within a viewer'''
      def __init__(self, viewer, accesscmd):
        self.accesscmd = accesscmd
        self.viewer = viewer
        
      def __getattr__(self,name):
        '''auto-instantiate javascript calls through model'''
        if name.startswith('_'): #object to ipython canary functions
            raise AttributeError("%r object has no attribute %r" %  (self.__class__, name))

        def makejs(*args,**kwargs):
          cmd = '\t%s.%s(' % (self.accesscmd,name);
          for arg in args:
              cmd += '%s,' % tostr(arg)
          cmd = cmd.rstrip(',')
          cmd += ');\n';

          self.viewer.startjs += cmd
          self.viewer.updatejs += cmd
          return self.viewer
        return makejs #return from getattr

    def getModel(self,*args,**kwargs):
      if self.viewergrid:
        if kwargs and 'viewer' in kwargs:
          coords = kwargs['viewer']
          if len(coords) != 2 or coords[0] >= self.viewergrid[0] or coords[1] >= self.viewergrid[1] or coords[0] < 0 or coords[1] < 0:
              raise ValueError("Incorrectly formated viewer argument.  Must specify row and column",self.viewergrid)
          cmd = '\tviewergrid_UNIQUEID[%d][%d].getModel(' % (coords[0],coords[1]);
          for arg in args:
              cmd += '%s,' % tostr(arg)
          cmd = cmd.rstrip(',')
          cmd += ')';
        else:
          raise ValueError('Must specify specific viewer with getModel and viewergrid enabled')
      else:
        cmd = '\tviewer_UNIQUEID.getModel(';
        for arg in args:
            cmd += '%s,' % tostr(arg)
        cmd = cmd.rstrip(',')
        cmd += ')';
      return self.model(self,cmd)

    def __getattr__(self,name):
        '''auto-instantiate javascript calls based on whatever the user provided'''
        if name.startswith('_'): #object to ipython canary functions
            raise AttributeError("%r object has no attribute %r" %  (self.__class__, name))
            
        def makejs(*args,**kwargs):
            if self.viewergrid:
                if kwargs and 'viewer' in kwargs:
                    coords = kwargs['viewer']
                    if len(coords) != 2 or coords[0] >= self.viewergrid[0] or coords[1] >= self.viewergrid[1] or coords[0] < 0 or coords[1] < 0:
                        raise ValueError("Incorrectly formated viewer argument.  Must specify row and column",self.viewergrid)
                    cmd = '\tviewergrid_UNIQUEID[%d][%d].%s(' % (coords[0],coords[1],name);
                    for arg in args:
                        cmd += '%s,' % tostr(arg)
                    cmd = cmd.rstrip(',')
                    cmd += ');\n';
                else: #apply to every viewer
                    cmd = ''
                    for r in range(self.viewergrid[0]):
                        for c in range(self.viewergrid[1]):
                            cmd += '\tviewergrid_UNIQUEID[%d][%d].%s(' % (r,c,name);
                            for arg in args:
                                cmd += '%s,' % tostr(arg)
                            cmd = cmd.rstrip(',')
                            cmd += ');\n';
            else:
                cmd = '\tviewer_UNIQUEID.%s(' % name;
                for arg in args:
                    cmd += '%s,' % tostr(arg)
                cmd = cmd.rstrip(',')
                cmd += ');\n';

            self.startjs += cmd
            self.updatejs += cmd
            return self

        return makejs
