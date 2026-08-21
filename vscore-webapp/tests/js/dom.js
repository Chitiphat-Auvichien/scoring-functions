var window = this; var self = this;
var navigator = { userAgent: "jsc" };
function El(){ this.style={}; this.children=[]; this.classList={add:function(){},remove:function(){}};
  this.appendChild=function(c){this.children.push(c);return c};
  this.addEventListener=function(){}; this.getContext=function(){return null};
  this.setAttribute=function(){}; this.getBoundingClientRect=function(){return {top:0,left:0,width:400,height:300}};
  this.querySelectorAll=function(){return []};
}
var document = { createElement:function(){return new El()}, createElementNS:function(){return new El()},
  body:new El(), addEventListener:function(){}, getElementById:function(){return new El()},
  querySelectorAll:function(){return []}, head:new El() };
function TextEncoder(){ this.encode=function(s){ var a=[]; for(var i=0;i<s.length;i++) a.push(s.charCodeAt(i)&255); return new Uint8Array(a); }; }
function TextDecoder(){ this.decode=function(b){ var s=""; for(var i=0;i<b.length;i++) s+=String.fromCharCode(b[i]); return s; }; }
function XMLHttpRequest(){ this.open=function(){}; this.send=function(){}; this.addEventListener=function(){}; }
var location = { href:"http://localhost/", protocol:"http:" };
var console = { log:function(){}, warn:function(){}, error:function(){}, info:function(){}, debug:function(){} };

var requestAnimationFrame = function(f){ f(); return 1; };
