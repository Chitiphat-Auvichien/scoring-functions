// Minimal DOM stub sufficient to run app.js and observe what it does.
var CALLS = { vibrate: null, addAtoms: null, setStyle: null, animate: null };
function El(id, tag) {
  this.id=id; this.tagName=(tag||"DIV").toUpperCase(); this.children=[];
  this.classList={_s:{},add:function(c){this._s[c]=1},remove:function(c){delete this._s[c]},
                  toggle:function(c,v){ v?this._s[c]=1:delete this._s[c] },
                  contains:function(c){return !!this._s[c]}};
  this.dataset={}; this.style={}; this.checked=true; this.value=""; this.hidden=false;
  this.textContent=""; this._listeners={};
  // innerHTML = "" must actually empty the node, as it does in a real DOM --
  // otherwise render()'s clear is a no-op and row counts accumulate.
  this._html="";
  Object.defineProperty(this,"innerHTML",{
    get:function(){return this._html;},
    set:function(v){ this._html=v; if(v==="") this.children.length=0; }
  });
  this.addEventListener=function(k,f){(this._listeners[k]=this._listeners[k]||[]).push(f);};
  this.fire=function(k,e){ (this._listeners[k]||[]).forEach(function(f){f(e||{});}); };
  this.appendChild=function(c){this.children.push(c); return c};
  this._qcache={};
  this.querySelector=function(sel){
    if(/pick/.test(sel)){ if(!this._pick){ this._pick=new El("input"); this._pick.className="pick";
      this._pick.dataset={i:String(this.dataset.index)}; } return this._pick; }
    // stable per selector: a fresh node each call means a listener attached
    // here is never the node a later lookup returns
    if(!this._qcache[sel]) this._qcache[sel]=new El("q");
    return this._qcache[sel]; };
  this.querySelectorAll=function(){return []};
  this.scrollIntoView=function(){};
  this.remove=function(){};
  this.click=function(){ this.fire('click',{target:{classList:{contains:function(){return false}}}}); };
}
var REG={};
var ABSENT={};
function get(id){ if(ABSENT[id]) return null; if(!REG[id]) REG[id]=new El(id); return REG[id]; }

REG["payload"]=new El("payload"); REG["payload"].textContent=PAYLOAD_TXT;
REG["downloads"]=new El("downloads"); REG["downloads"].textContent=DL_TXT;
var TBODY=new El("tbody","TBODY");
var document={
  getElementById:get,
  createElement:function(t){return new El("new",t)},
  querySelector:function(s){ return s.indexOf("tbody")>-1 ? TBODY : new El("s"); },
  querySelectorAll:function(s){
    if (s.indexOf("tbody tr")>-1) return TBODY.children;
    return [];
  },
  addEventListener:function(){},
  body:new El("body")
};
var window={matchMedia:function(){return{matches:false}}};
var URL={createObjectURL:function(){return "blob:x"},revokeObjectURL:function(){}};
function Blob(a,b){ this.parts=a; this.type=b.type; }
function getComputedStyle(){ return {getPropertyValue:function(){return "white"}}; }
// defaults the real inputs would carry
get("f-lo").value="0"; get("f-hi").value="99999";
get("amp").value="1"; get("frm").value="10";

get("bonds").checked=true;
get("axes").checked=true;

["compare","cmp-n","cmp-scores","cmp-bonds","cmp-clear"].forEach(function(id){ get(id); });

ABSENT["f-ref"]=true; delete REG["f-ref"];
