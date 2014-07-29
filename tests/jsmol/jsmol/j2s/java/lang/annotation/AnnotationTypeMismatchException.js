Clazz.load(["java.lang.RuntimeException"],"java.lang.annotation.AnnotationTypeMismatchException",null,function(){
c$=Clazz.decorateAsClass(function(){
this.$element=null;
this.$foundType=null;
Clazz.instantialize(this,arguments);
},java.lang.annotation,"AnnotationTypeMismatchException",RuntimeException);
Clazz.makeConstructor(c$,
function(element,foundType){
Clazz.superConstructor(this,java.lang.annotation.AnnotationTypeMismatchException,[("annotation.1",element,foundType)]);
this.$element=element;
this.$foundType=foundType;
},"java.lang.reflect.Method,~S");
Clazz.defineMethod(c$,"element",
function(){
return this.$element;
});
Clazz.defineMethod(c$,"foundType",
function(){
return this.$foundType;
});
});
