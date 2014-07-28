Clazz.load(["java.lang.RuntimeException"],"java.lang.annotation.IncompleteAnnotationException",null,function(){
c$=Clazz.decorateAsClass(function(){
this.$annotationType=null;
this.$elementName=null;
Clazz.instantialize(this,arguments);
},java.lang.annotation,"IncompleteAnnotationException",RuntimeException);
Clazz.makeConstructor(c$,
function(annotationType,elementName){
Clazz.superConstructor(this,java.lang.annotation.IncompleteAnnotationException,[("annotation.0",elementName,annotationType)]);
this.$annotationType=annotationType;
this.$elementName=elementName;
},"Class,~S");
Clazz.defineMethod(c$,"annotationType",
function(){
return this.$annotationType;
});
Clazz.defineMethod(c$,"elementName",
function(){
return this.$elementName;
});
});
