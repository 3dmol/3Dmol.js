import time
import json
import sys
import re
from urllib.parse import urlsplit, urlunsplit

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
    '''Construct HTML/JavaScript views for IPython notebooks and exported pages.
       Displayed views remain interactive without a running kernel. Python queues
       browser commands; JavaScript results and events are not returned to Python.
       The page still needs JavaScript, WebGL, and any referenced network resources.

       A viewergrid tuple (rows, columns) creates a grid sharing one canvas.
       Commands apply to all cells unless viewer=(row, column) selects one cell
       using zero-based indices. getModel requires a cell when a grid is enabled.

       Forwarded methods use the GLViewer API with JSON-serializable positional
       arguments and return this Python view for chaining, not JavaScript results.
       getModel returns a proxy for queuing GLModel calls. Python callbacks and
       browser-side promises are not supported across this one-way interface.
       https://3dmol.org/doc/GLViewer.html
    '''
    def __init__(self,query='',width=640,height=480,viewergrid=None,data=None,style=None,linked=True,options=dict(),format=None,js=None,
                 edit=False,trajectory=False,learning=False,ui_js=None):
        '''Create a 3Dmol.js view.
            width -- container width in pixels (int) or a CSS dimension string
            height -- container height in pixels (int) or a CSS dimension string
            query -- optional argument to provide to $3Dmol.download
            viewergrid -- optional tuple (rows,columns) to define grid
            linked -- whether mouse controls affect all grid cells (default True)
            data -- molecular data to provide to addModel, with viewer grid can be indexed (r,c)
            format -- format of provided data
            style -- style to apply, with viewer grid can be indexed (r,c)
            options -- optional options to provide to $3Dmol.download or addModel
            js -- core bundle URL; defaults to the same pinned npm release for all
                  views, with or without controls. The first load wins.
            edit -- show the edit-menu button (default False)
            trajectory -- show trajectory controls when frames are present (default False)
            learning -- True for a classroom button, or a LearningConfig dictionary;
                        hidden mode requires {'mode': 'hidden', 'session': 'unique-id'}.
                        The server defaults to https://3dmol.csb.pitt.edu.
            ui_js -- optional UI bundle URL, otherwise inferred alongside js;
                     use matching core and UI builds when overriding either URL
            Controls require a single viewer (viewergrid=(1,1) is also supported).
            Python queues commands; classroom events do not call back into Python.'''
        if not isinstance(learning, (bool, dict)):
            raise TypeError('learning must be a boolean or a configuration dictionary')
        if isinstance(learning, dict):
            mode = learning.get('mode', 'button')
            if mode not in ('button', 'hidden'):
                raise ValueError("learning mode must be 'button' or 'hidden'")
            session = learning.get('session')
            if mode == 'hidden' and (not isinstance(session, str) or not session.strip()):
                raise ValueError('Hidden classroom mode requires a nonblank session identifier')
        has_ui = bool(edit or trajectory or learning is True or isinstance(learning, dict))
        if has_ui and viewergrid and tuple(viewergrid) != (1, 1):
            raise ValueError('Controls require a single viewer; use separate views instead of a multi-cell viewergrid')
        if js is None:
            js = 'https://cdn.jsdelivr.net/npm/3dmol@3.0.0/build/3Dmol-min.js'
        if has_ui and ui_js is None:
            parts = urlsplit(js)
            path, matches = re.subn(r'(^|/)3Dmol(-min)?\.js$', r'\g<1>3Dmol.ui\g<2>.js', parts.path)
            if not matches:
                raise ValueError('Specify ui_js when js does not end in 3Dmol.js or 3Dmol-min.js')
            ui_js = urlunsplit(parts._replace(path=path))
        # These values live inside an inline script, not just a JavaScript string.
        # Escape '<' so a session identifier or URL cannot close the script tag.
        ui_config = json.dumps({'edit': bool(edit), 'trajectory': bool(trajectory),
                                'learning': learning}).replace('<', '\\u003c')
        divid = "3dmolviewer_UNIQUEID"
        warnid = "3dmolwarning_UNIQUEID"
        self.uniqueid = None
        #convert numerical width/height to pixels
        if type(width) == int:
            width = '%dpx'%width
        if type(height) == int:
            height = '%dpx'%height
        self.startjs = '''<div id="%s"  style="position: relative; width: %s; height: %s;">
        <p id="%s" style="background-color:#ffcccc;color:black">3Dmol.js failed to load for some reason.  Please check your browser console for error messages.<br></p>
        </div>\n''' % (divid,width,height,warnid)
        self.startjs += '<script>\n'
        self.endjs = '</script>'

        self.updatejs = '' # code added since last show
        #load 3dmol, but only once, can't use jquery :-(
        #https://medium.com/@vschroeder/javascript-how-to-execute-code-from-an-asynchronously-loaded-script-although-when-it-is-not-bebcbd6da5ea
        self.startjs += """
var loadScriptAsync = function(uri){
  return new Promise((resolve, reject) => {
    //this is to ignore the existence of requirejs amd
    var savedexports, savedmodule;
    if (typeof exports !== 'undefined') savedexports = exports;
    else exports = {}
    if (typeof module !== 'undefined') savedmodule = module;
    else module = {}

    var tag = document.createElement('script');
    tag.src = uri;
    tag.async = true;
    tag.onload = () => {
        exports = savedexports;
        module = savedmodule;
        resolve();
    };
    tag.onerror = () => {
        exports = savedexports;
        module = savedmodule;
        reject(new Error('Unable to load ' + uri));
    };
  var firstScriptTag = document.getElementsByTagName('script')[0];
  firstScriptTag.parentNode.insertBefore(tag, firstScriptTag);
});
};

if(typeof $3Dmolpromise === 'undefined') {
  $3Dmolpromise = typeof $3Dmol === 'undefined' ? loadScriptAsync(%s) : Promise.resolve();
}

var viewer_UNIQUEID = null;
""" % json.dumps(js).replace('<', '\\u003c')

        dependency = '$3Dmolpromise'
        if has_ui:
            self.startjs += """
if(typeof $3DmolUIpromise === 'undefined') {
  $3DmolUIpromise = $3Dmolpromise.then(function() {
    if (typeof $3Dmol.StateManager === 'function') return;
    return loadScriptAsync(%s);
  }).then(function() {
    if (typeof $3Dmol.StateManager !== 'function') {
      throw new Error('The UI bundle does not provide StateManager. Use matching core and UI builds.');
    }
  });
}
""" % json.dumps(ui_js).replace('<', '\\u003c')
            dependency = '$3DmolUIpromise'

        self.startjs += 'var viewer_ready_UNIQUEID = %s.then(function() {\n' % dependency
        self.startjs += """
var warn = document.getElementById("%s");
if(warn) {
    warn.parentNode.removeChild(warn);
}
""" % warnid

        self.endjs = '''});
viewer_ready_UNIQUEID.catch(function(error) {
    var host = document.getElementById("%s");
    if (host) {
        var warning = document.getElementById("%s") || document.createElement('p');
        warning.id = "%s";
        warning.style.cssText = 'position:absolute;top:0;background:#ffcccc;color:black';
        warning.textContent = '3Dmol.js failed to initialize: ' + error.message;
        host.appendChild(warning);
    }
    console.error(error);
});
</script>''' % (divid, warnid, warnid)
        self.viewergrid = None

        if viewergrid: #create split view
            if len(viewergrid) != 2:
                raise ValueError("Incorrectly formated viewergrid arguments.  Must specify rows x columns",viewergrid)
            self.startjs += "var viewergrid_UNIQUEID = null;\n";
            self.startjs += 'viewergrid_UNIQUEID = $3Dmol.createViewerGrid(document.getElementById("%s"),{rows: %d, cols: %d, control_all: %s},{backgroundColor:"white"});\n' % (divid, viewergrid[0],viewergrid[1],'true' if linked else 'false')
            self.startjs += "viewer_UNIQUEID = viewergrid_UNIQUEID[0][0];\n"
            self.viewergrid = viewergrid
        else:
            self.startjs += 'viewer_UNIQUEID = $3Dmol.createViewer(document.getElementById("%s"),{backgroundColor:"white"});\n' % divid
        if query:
            if viewergrid:
                for r in range(viewergrid[0]):
                    for c in range(viewergrid[1]):
                        self.startjs += 'return $3Dmol.download(%s, viewergrid_UNIQUEID[%d][%d], %s).then(function() {\n' % (json.dumps(query),r,c,json.dumps(options))
                        self.endjs = "});\n" + self.endjs
            else:
                self.startjs += 'return $3Dmol.download(%s, viewer_UNIQUEID, %s).then(function() {\n' % (json.dumps(query),json.dumps(options))
                self.endjs = "});\n" + self.endjs

        finish = ''
        if has_ui:
            if learning is True or isinstance(learning, dict):
                finish += "if (typeof $3Dmol.LearningSession !== 'function') throw new Error('This UI bundle does not support classrooms. Use matching current core and UI builds.');\n"
            finish += 'new $3Dmol.StateManager(viewer_UNIQUEID, %s);\n' % ui_config
        self.endjs = finish + 'return viewer_UNIQUEID;\n' + self.endjs
                
        if viewergrid:
            for r in range(viewergrid[0]):
                for c in range(viewergrid[1]):
                    cmds = ''
                    if data:
                        try:
                            if type(data) == str:
                              d = data
                            else:
                              d = data[r][c]
                        except:
                            d = data
                        try:
                            if type(format) == str:
                              f = format
                            else:
                              f = format[r][c]
                        except:
                            f = format
                        self.startjs += f"viewergrid_UNIQUEID[{r}][{c}].addModel({json.dumps(d)},{json.dumps(f) if f else 'undefined'},{json.dumps(options) if options else 'undefined'});\n"
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
                cmds = f"viewer_UNIQUEID.addModel({json.dumps(data)},{json.dumps(format) if format else 'undefined'},{json.dumps(options) if options else 'undefined'});\n"                
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
        html += '''<script>document.getElementById("%s").append(document.getElementById("3dmolviewer_%s")); </script>'''%(containerid,self.uniqueid)
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
            viewer_ready_%s.then(function() { //wait for this view's scene and controls
                %s
                viewer_%s.render();
                if (viewer_%s.ui) viewer_%s.ui.refresh();
            });
            </script>''' % (self.uniqueid,self.updatejs.replace('UNIQUEID',self.uniqueid),self.uniqueid,self.uniqueid,self.uniqueid)
        self.updatejs = ''
        return IPython.display.publish_display_data({'application/3dmoljs_load.v0':script, 'text/html': script},metadata={})


    def write_html(self, f=None, fullpage=False):
      '''Export the queued scene as HTML without requiring notebook display.
      f -- file name (str) or writable file object; omitted returns an HTML fragment
      fullpage -- fill the page instead of using the configured size when writing f;
                  this option is not applied when returning a string
      Referenced scripts and remote molecular data are not embedded by this method.
      '''
      
      if f == None:
        return self._make_html()
        
      if type(f) == str:
        f = open(f,'wt')
      html = self._make_html()
      html = f'''<html>
<body style="margin: 0; padding: 0; display: block;">
{html}
</body>
</html>'''
      if fullpage:
        html = re.sub(r'width: (\S+);', 'width: 100%;', html)
        html = re.sub(r'height: (\S+);', 'height: 100vh;', html)
        
      f.write(html)
      
      
    @using_ipython
    def png(self):
        '''Display a PNG of an already instantiated viewer in notebook output.
        The image is generated in the browser; no image bytes are returned to Python.
        '''
        if not self.uniqueid:
            raise AssertionError('Must instantiate viewer before generating image.')
        script = '''<img id="img_{0}">
            <script>
            viewer_ready_{0}.then(function() {{
            var png = viewer_{0}.pngURI()
            document.getElementById('img_{0}').src = png;
            }});
            </script>'''.format(self.uniqueid)
        return IPython.display.publish_display_data({'application/3dmoljs_load.v0':script, 'text/html': script},metadata={})
    
    @using_ipython
    def apng(self, nframes=1):
        '''Display an APNG of an already instantiated viewer in notebook output.
        Capture nframes rendered changes; no image bytes are returned to Python.
        '''
        if not self.uniqueid:
            raise AssertionError('Must instantiate viewer before generating animated image.')
        script = '''<img id="img_{0}">
            <script>
            viewer_ready_{0}.then(function() {{
            viewer_{0}.apngURI({1}).then(png => {{
            document.getElementById('img_{0}').src = png; }});
            }});
            </script>'''.format(self.uniqueid,nframes)
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
      """Return a proxy for queuing GLModel methods, not its atom data.
      Positional arguments are forwarded to GLViewer.getModel. With a viewer grid,
      viewer=(row, column) is required. Proxy calls return the parent Python view.
      """
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
        if name.startswith('_') or name == 'getdoc': #object to ipython canary functions
            raise AttributeError("%r object has no attribute %r" %  (self.__class__, name))
            
        print_result = False
        if name.startswith('print_to_console_'):
            print_result = True
            name = name[17:]

        def makejs(*args,**kwargs):

            def make_viewer_cmd(vname, name, *args, **kwargs):
                cmd = '\t'
                if print_result:
                    cmd += 'console.log('
                cmd += f'{vname}.{name}('
                for arg in args:
                    cmd += '%s,' % tostr(arg)
                cmd = cmd.rstrip(',')
                if print_result:
                    cmd += ')'                         
                cmd += ');\n'
       
                return cmd

            if self.viewergrid:
                if kwargs and 'viewer' in kwargs:
                    coords = kwargs['viewer']
                    if len(coords) != 2 or coords[0] >= self.viewergrid[0] or coords[1] >= self.viewergrid[1] or coords[0] < 0 or coords[1] < 0:
                        raise ValueError("Incorrectly formated viewer argument.  Must specify row and column",self.viewergrid)
                    cmd = make_viewer_cmd('viewergrid_UNIQUEID[%d][%d]'%(coords[0],coords[1]),name,*args,**kwargs)
                else: #apply to every viewer
                    cmd = ''
                    for r in range(self.viewergrid[0]):
                        for c in range(self.viewergrid[1]):
                            cmd += make_viewer_cmd('viewergrid_UNIQUEID[%d][%d]'%(r,c),name,*args,**kwargs)
            else:
                cmd = make_viewer_cmd('viewer_UNIQUEID',name,*args,**kwargs)

            self.startjs += cmd
            self.updatejs += cmd
            if print_result:
                self.update()
                return "Inspect the JavaScript console for your requested result."
            else:
                return self

        return makejs
