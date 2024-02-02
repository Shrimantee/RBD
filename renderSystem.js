

o3djs.base.o3d = o3d;
o3djs.require('o3djs.webgl');
o3djs.require('o3djs.math');
o3djs.require('o3djs.rendergraph');
o3djs.require('o3djs.primitives');
o3djs.require('o3djs.material');

var g_o3dElement;
var g_client;
var g_o3d;
var g_math;
var g_pack;
var g_viewInfo;
var g_thisRenderer;

Renderer = function ()
{}

Renderer.prototype.initialize = function (eyePosition, targetPosition, displayOutput, bodyColors)
{
    this.eyePosition = eyePosition;
    this.targetPosition = targetPosition;
    this.displayOutput = displayOutput;
    this.bodyTransforms = [];
    this.markerScale = 0.1;
    this.bodyColors = bodyColors;
    g_thisRenderer = this;
    o3djs.webgl.makeClients(Renderer.prototype.initialize2);
}

Renderer.prototype.initialize2 = function (clientElements)
{
    g_thisRenderer.initGlobals(clientElements);
    g_thisRenderer.initContext();
    g_thisRenderer.createShapes();
    g_client.setRenderCallback(this.displayOutput);
}

Renderer.prototype.initGlobals = function (clientElements)
{
    g_o3dElement = clientElements[0];
    window.g_client = g_client = g_o3dElement.client;
    g_o3d = g_o3dElement.o3d;
    g_math = o3djs.math;
    
    g_pack = g_client.createPack();
    
    var clearColor = [.98, .98, .98, 1];
    g_viewInfo = o3djs.rendergraph.createBasicView(
        g_pack, g_client.root, g_client.renderGraphRoot, clearColor);
}


Renderer.prototype.initContext = function ()
{
    g_viewInfo.drawContext.projection = g_math.matrix4.perspective(
        g_math.degToRad(50), g_o3dElement.clientWidth / g_o3dElement.clientHeight, 1, 5000);
    g_viewInfo.drawContext.view = g_math.matrix4.lookAt(
        this.eyePosition, this.targetPosition, [0, 1, 0]);
}

Renderer.prototype.createShapes = function ()
{
    
    
    
    for (var b in g_world.rigidBodies)
    {

        var color = (this.bodyColors) ? this.bodyColors[b] : [1.0, 0.2, 1.0, 1];
        
        var materialCube = o3djs.material.createBasicMaterial(
            g_pack, g_viewInfo, color);
        materialCube.getParam('specularFactor').value = 0.2;
        materialCube.getParam('ambient').value = [Math.random(), Math.random(),Math.random(), 1];
        console.log(""+b)
        // if(b>8 || b==11){
            var cube = o3djs.primitives.createCube(
                g_pack, materialCube, 1.0);
        // }
        // else if(b<8 ){
        //     var cube = o3djs.primitives.createCylinder(
        //         g_pack, materialCube, 0.5,1,20,20);
        // }
        
        var transform = g_pack.createObject('Transform');
        transform.addShape(cube);
        transform.parent = g_client.root;
        this.bodyTransforms.push(transform);
    }

    
}


Renderer.prototype.update = function ()
{
    for (var bi in g_world.rigidBodies)
    {
        var b = g_world.rigidBodies[bi];
        var t = this.bodyTransforms[bi];
        t.identity();
        t.translate(b.x.x, b.x.y, b.x.z);
        o3d.Transform.compose(t.localMatrix, b.R.getAsO3DMatrix4());
        t.scale(b.l.x, b.l.y, b.l.z);
    }
   
}