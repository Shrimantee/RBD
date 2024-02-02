
/* UTILITY METHODS */
function Vector(x, y, z)
{
    this.x = x;
    this.y = y;
    this.z = z;
}
Vector.prototype.size = function ()
{
    return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
}

Vector.prototype.negate = function ()
{
    this.x = -this.x;
    this.y = -this.y;
    this.z = -this.z;
}
Vector.prototype.add = function (v)
{
    return new Vector(this.x + v.x, this.y + v.y, this.z + v.z);
}

Vector.prototype.sub = function (v)
{
    return new Vector(this.x - v.x, this.y - v.y, this.z - v.z);
}




Vector.prototype.scalarMult = function (a)
{
    return new Vector(this.x * a, this.y * a, this.z * a);
}
Vector.prototype.vecMult = function (v)
{
    return new Vector(this.x * v.x, this.y * v.y, this.z * v.z);
}
Vector.prototype.scalarDiv = function (a)
{
    return new Vector(this.x / a, this.y / a, this.z / a);
}

Vector.prototype.cross = function (v)
{
    return new Vector(this.y * v.z - this.z * v.y, 
                      this.z * v.x - this.x * v.z,
                      this.x * v.y - this.y * v.x);
}

Vector.prototype.accumulate = function (v)
{
    this.x += v.x;
    this.y += v.y;
    this.z += v.z;
}

Vector.prototype.degeneratw = function (v)
{
    this.x -= v.x;
    this.y -= v.y;
    this.z -= v.z;
}
Vector.prototype.dot = function (v)
{
    return this.x * v.x + this.y * v.y + this.z * v.z;
}

Vector.prototype.copy = function (a)
{
    return new Vector(a.x, a.y, a.z);
}


Vector.prototype.norm = function ()
{
    var l = Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
    this.x /= l;
    this.y /= l;
    this.z /= l;
}



Vector.prototype.getSkewSymmMatrix = function ()
{
    var m = new Matrix2DH();
    m.coeff00 = 0;
    m.coeff01 = -this.z;
    m.coeff02 = this.y;
    m.coeff10 = this.z;
    m.coeff11 = 0;
    m.coeff12 = -this.x;
    m.coeff20 = -this.y;
    m.coeff21 = this.x;
    m.coeff22 = 0;
    return m;
}

var Vec = new Vector(0,0,0);

function Matrix2DH()
{}

Matrix2DH.prototype.setIdMat = function ()
{
    this.coeff00 = 1;this.coeff01 = 0;this.coeff02 = 0;
    this.coeff10 = 0;this.coeff11 = 1;this.coeff12 = 0;
    this.coeff20 = 0;this.coeff21 = 0;this.coeff22 = 1;
}

Matrix2DH.prototype.transpose = function ()
{
    var newMat = new Matrix2DH();
    newMat.coeff00 = this.coeff00;
    newMat.coeff01 = this.coeff10;
    newMat.coeff02 = this.coeff20;
    newMat.coeff10 = this.coeff01;
    newMat.coeff11 = this.coeff11;
    newMat.coeff12 = this.coeff21;
    newMat.coeff20 = this.coeff02;
    newMat.coeff21 = this.coeff12;
    newMat.coeff22 = this.coeff22;
    return newMat;
}

Matrix2DH.prototype.copy = function (mat)
{
    var newMat = new Matrix2DH();
    newMat.coeff00 = mat.coeff00;newMat.coeff01 = mat.coeff01;newMat.coeff02 = mat.coeff02;
    newMat.coeff10 = mat.coeff10;newMat.coeff11 = mat.coeff11;newMat.coeff12 = mat.coeff12;
    newMat.coeff20 = mat.coeff20;newMat.coeff21 = mat.coeff21;newMat.coeff22 = mat.coeff22;
    return newMat;
}



Matrix2DH.prototype.accumulate = function (mat)
{
    this.coeff00 += mat.coeff00;
    this.coeff01 += mat.coeff01;
    this.coeff02 += mat.coeff02;
    this.coeff10 += mat.coeff10;
    this.coeff11 += mat.coeff11;
    this.coeff12 += mat.coeff12;
    this.coeff20 += mat.coeff20;
    this.coeff21 += mat.coeff21;
    this.coeff22 += mat.coeff22;
}

Matrix2DH.prototype.mul = function (mat)
{
    var newMat = new Matrix2DH();
    newMat.coeff00 = this.coeff00 * mat.coeff00 + this.coeff01 * mat.coeff10 + this.coeff02 * mat.coeff20;
    newMat.coeff01 = this.coeff00 * mat.coeff01 + this.coeff01 * mat.coeff11 + this.coeff02 * mat.coeff21;
    newMat.coeff02 = this.coeff00 * mat.coeff02 + this.coeff01 * mat.coeff12 + this.coeff02 * mat.coeff22;

    newMat.coeff10 = this.coeff10 * mat.coeff00 + this.coeff11 * mat.coeff10 + this.coeff12 * mat.coeff20;
    newMat.coeff11 = this.coeff10 * mat.coeff01 + this.coeff11 * mat.coeff11 + this.coeff12 * mat.coeff21;
    newMat.coeff12 = this.coeff10 * mat.coeff02 + this.coeff11 * mat.coeff12 + this.coeff12 * mat.coeff22;

    newMat.coeff20 = this.coeff20 * mat.coeff00 + this.coeff21 * mat.coeff10 + this.coeff22 * mat.coeff20;
    newMat.coeff21 = this.coeff20 * mat.coeff01 + this.coeff21 * mat.coeff11 + this.coeff22 * mat.coeff21;
    newMat.coeff22 = this.coeff20 * mat.coeff02 + this.coeff21 * mat.coeff12 + this.coeff22 * mat.coeff22;

    return newMat;
}

Matrix2DH.prototype.vecMult = function (vec)
{
    return new Vector(this.coeff00 * vec.x + this.coeff01 * vec.y + this.coeff02 * vec.z,
                      this.coeff10 * vec.x + this.coeff11 * vec.y + this.coeff12 * vec.z,
                      this.coeff20 * vec.x + this.coeff21 * vec.y + this.coeff22 * vec.z);
}

Matrix2DH.prototype.scalarMult = function (scale)
{
    var newMat = new Matrix2DH();
    newMat.coeff00 = this.coeff00 * scale;
    newMat.coeff01 = this.coeff01 * scale;
    newMat.coeff02 = this.coeff02 * scale;
    newMat.coeff10 = this.coeff10 * scale;
    newMat.coeff11 = this.coeff11 * scale;
    newMat.coeff12 = this.coeff12 * scale;
    newMat.coeff20 = this.coeff20 * scale;
    newMat.coeff21 = this.coeff21 * scale;
    newMat.coeff22 = this.coeff22 * scale;
    return newMat;
}
Matrix2DH.prototype.add = function (mat)
{
    var newMat = new Matrix2DH();
    newMat.coeff00 = this.coeff00 + mat.coeff00;
    newMat.coeff01 = this.coeff01 + mat.coeff01;
    newMat.coeff02 = this.coeff02 + mat.coeff02;
    newMat.coeff10 = this.coeff10 + mat.coeff10;
    newMat.coeff11 = this.coeff11 + mat.coeff11;
    newMat.coeff12 = this.coeff12 + mat.coeff12;
    newMat.coeff20 = this.coeff20 + mat.coeff20;
    newMat.coeff21 = this.coeff21 + mat.coeff21;
    newMat.coeff22 = this.coeff22 + mat.coeff22;
    return newMat;
}
Matrix2DH.prototype.getAsArray4 = function ()
{
    return [this.coeff00, this.coeff10, this.coeff20, 0, 
            this.coeff01, this.coeff11, this.coeff21, 0, 
            this.coeff02, this.coeff12, this.coeff22, 0,
            0, 0, 0, 1];
}

Matrix2DH.prototype.getAsO3DMatrix4 = function ()
{
    return [[this.coeff00, this.coeff10, this.coeff20, 0], 
            [this.coeff01, this.coeff11, this.coeff21, 0], 
            [this.coeff02, this.coeff12, this.coeff22, 0], 
            [0, 0, 0, 1]];
}

var mat3D = new Matrix2DH();


function Body()
{

    this.force = new Vector(0, 0, 0);
    this.torque = new Vector(0, 0, 0);
}

Body.prototype.set = function (rb)
{
    this.mass     = rb.mass;
    this.Ibody    = mat3D.copy(rb.Ibody);
    this.Ibodyinv = mat3D.copy(rb.Ibodyinv);
    this.l        = Vec.copy(rb.l);

    this.x        = Vec.copy(rb.x);
    this.R        = mat3D.copy(rb.R);
    this.P        = Vec.copy(rb.P);
    this.L        = Vec.copy(rb.L);
    this.Iinv     = mat3D.copy(rb.Iinv);
    this.v        = Vec.copy(rb.v);
    this.omega    = Vec.copy(rb.omega);
    this.force    = Vec.copy(rb.force);
    this.torque   = Vec.copy(rb.torque);
}

Body.prototype.applyAuxi = function ()
{
    this.v = this.P.scalarDiv(this.mass);
    this.Iinv = this.R.mul(this.Ibodyinv).mul(this.R.transpose());
    this.omega = this.Iinv.vecMult(this.L);
}

Body.prototype.integrate = function (dt)
{
/*     this.x.accumulate(this.v.scalarMult(dt));
    var m = this.omega.getSkewSymmMatrix().mul(this.R);
    this.R.accumulate(m.scalarMult(dt));

    this.P.accumulate(this.force.scalarMult(dt));
    this.L.accumulate(this.torque.scalarMult(dt));

    var kdf = this.v.scalarMult(g_world.kdl * dt);
    kdf.negate();
    this.P.accumulate(kdf);
    var kdt = this.omega.scalarMult(g_world.kdw * dt);
    kdt.negate();
    this.L.accumulate(kdt); */
    this.x.accumulate(this.v.scalarMult(dt));
    var m = this.omega.getSkewSymmMatrix().mul(this.R);
    this.R.accumulate(m.scalarMult(dt));

    this.P.accumulate(this.force.scalarMult(dt));
    this.L.accumulate(this.torque.scalarMult(dt));

    var kdf = this.v.scalarMult(g_world.kdl * dt);
    kdf.negate();
    this.P.accumulate(kdf);
    var kdt = this.omega.scalarMult(g_world.kdw * dt);
    kdt.negate();
    this.L.accumulate(kdt);
}

Body.prototype.renormalizeRotationMatrix = function ()
{
    var v0 = new Vector(this.R.coeff00, this.R.coeff10, this.R.coeff20);
    v0.norm();
    var v1 = new Vector(this.R.coeff01, this.R.coeff11, this.R.coeff21);
    v1.norm();
    var v2 = v0.cross(v1);
    v1 = v2.cross(v0);
    this.R.coeff00 = v0.x;
    this.R.coeff01 = v1.x;
    this.R.coeff02 = v2.x;
    this.R.coeff10 = v0.y;
    this.R.coeff11 = v1.y;
    this.R.coeff12 = v2.y;
    this.R.coeff20 = v0.z;
    this.R.coeff21 = v1.z;
    this.R.coeff22 = v2.z;
}

Body.prototype.createRB = function (lx, ly, lz, m, x, R)
{
    this.l = new Vector(lx, ly, lz);
    this.mass = m;
    this.x = x;
    this.R = R;
    this.P = new Vector(0, 0, 0);
    this.L = new Vector(0, 0, 0);
    this.Ibody = new Matrix2DH();
    this.Ibody.setIdMat();
    var m0 = this.mass / 12;
    // console.log(""+this.mass)

    this.Ibody.coeff00 = m0 * (ly * ly + lz * lz);
    this.Ibody.coeff11 = m0 * (lx * lx + lz * lz);
    this.Ibody.coeff22 = m0 * (lx * lx + ly * ly);
    this.Ibodyinv = new Matrix2DH();
    this.Ibodyinv.setIdMat();
    this.Ibodyinv.coeff00 = 1 / this.Ibody.coeff00;
    this.Ibodyinv.coeff11 = 1 / this.Ibody.coeff11;
    this.Ibodyinv.coeff22 = 1 / this.Ibody.coeff22;
    this.applyAuxi();
}

function PhysicsWorld()
{
    this.gravity = new Vector(0, -5.8, 0);
    this.rigidBodies = [];
    this.prevRB = [];
    this.steps = 20;
    this.kdl = 0.1; // linear 
    this.kdw = 0.1; // angular 
}

var g_world = new PhysicsWorld();

PhysicsWorld.prototype.ptVelocity = function (b, p)
{
    return b.v.add(b.omega.cross(p.sub(b.x))); 
}

PhysicsWorld.prototype.copyBodies = function (ba0, ba1, reuse) 
{
    if (reuse)
    {
        for (var i = 0; i < ba1.length; i++)
        {
            ba0[i].setFrom(ba1[i]);
        }
    }
    else
    {
        ba0.length = 0;
        for (var i = 0; i < ba1.length; i++)
        {
            var b0 = new Body();
            b0.set(ba1[i]);
            ba0.push(b0);
        }
    }
}
PhysicsWorld.prototype.integrateBodies = function (dt)
{
    for (var i = 0; i < this.rigidBodies.length; i++)
    {
        var body = this.rigidBodies[i];
        if (body.mass >= 0)
        {
            body.integrate(dt);
            body.renormalizeRotationMatrix();
            body.applyAuxi();
        }
    }
}

function Collision()
{
    this.depth = 0;
}

var boxPoints  = [ new Vector(-0.5, 0, 0), new Vector(0.5, 0, 0),
                    new Vector(0, -0.5, 0), new Vector(0, 0.5, 0),
                    new Vector(0, 0, -0.5), new Vector(0, 0, 0.5)];
var boxNorm = [ new Vector(-1, 0, 0),   new Vector(1, 0, 0),
                    new Vector(0, -1, 0),   new Vector(0, 1, 0),
                    new Vector(0, 0, -1),   new Vector(0, 0, 1)];

PhysicsWorld.prototype.linePlaneIntersection = function (p, n, pa, pb)
{
    var den = ((pb.sub(pa)).dot(n));
    if (den == 0)
        return 1.0e30; // manipulate to fine result
    return ((p.sub(pa)).dot(n))/den;
}

PhysicsWorld.prototype.detectCollisionBodyBody = function (bodyb, bodya0, bodya1)
{
    var contacts = [];
    var invRb = bodyb.R.transpose();
    var invlb = new Vector(1/bodyb.l.x, 1/bodyb.l.y, 1/bodyb.l.z);
    var reldist0 = bodya0.x.sub(bodyb.x);
    var reldist1 = bodya1.x.sub(bodyb.x);
    for (var x = -0.5; x <= 0.5; x += 1)
        for (var y = -0.5; y <= 0.5; y += 1)
            for (var z = -0.5; z <= 0.5; z += 1)
            {
                var pb = new Vector(x, y, z);
                pb = pb.vecMult(bodya1.l);                   
                var rpb = bodya1.R.vecMult(pb);          
                var wrpb = rpb.add(reldist1);             
                var qb = invRb.vecMult(wrpb);            
                qb = qb.vecMult(invlb);                 

                if ((qb.x < 0.5 && qb.x > -0.5) &&
                    (qb.y < 0.5 && qb.y > -0.5) &&
                    (qb.z < 0.5 && qb.z > -0.5))
                {
                    var con = new Collision();
                    //console.log("Collided");

                    var pa = new Vector(x, y, z);
                    pa = pa.vecMult(bodya0.l);             
                    var rpa = bodya0.R.vecMult(pa);      
                    var wrpa = rpa.add(reldist0);         
                    var qa = invRb.vecMult(wrpa);        
                    qa = qa.vecMult(invlb);            

                    var max_d = -1.0e30;
                    var max_f = -1;
                    for (var f = 0; f < 6; f++)
                    {
                        
                        var d = this.linePlaneIntersection(boxPoints[f], boxNorm[f], qa, qb);
                        if (d <= 1)
                        {
                            if (d > max_d)
                            {
                                max_d = d;
                                max_f = f;
                            }
                        }
                    }
                    var worldPos = qa.add((qb.sub(qa)).scalarMult(max_d));  
                    worldPos = worldPos.vecMult(bodyb.l);               
                    worldPos = bodyb.R.vecMult(worldPos);                
                    con.pos = worldPos.add(bodyb.x);                       
                    con.normal = boxNorm[max_f];
                    if(!qb)
                    console.log(""+qb)
                    con.depth = -((qb.sub(boxPoints[max_f])).dot(con.normal));
                    con.bodyb = bodyb;
                    con.bodya = bodya1;
                    contacts.push(con);
                }
            }
    return contacts;
}

PhysicsWorld.prototype.findCollisionPoints = function ()
{
    var contacts = [];
    for (var i = 0; i < this.rigidBodies.length; i++)
    {   /* console.log(this.rigidBodies.length) */
        var bodyb = this.rigidBodies[i];
        for (var j = 0; j < this.rigidBodies.length; j++)
        {
            if (i != j)
            {
                var bodya0 = this.prevRB[j];
                var bodya1 = this.rigidBodies[j];
                var c = this.detectCollisionBodyBody(bodyb, bodya0, bodya1);
                contacts.push.apply(contacts, c);
            }
        }
    }
    return contacts;
}

PhysicsWorld.prototype.handleCollision = function (con)
{
    var padot = this.ptVelocity(con.bodya, con.pos);
    var pbdot = this.ptVelocity(con.bodyb, con.pos);
    var ra = con.pos.sub(con.bodya.x);
    var rb = con.pos.sub(con.bodyb.x);
    var vrel = con.normal.dot(padot.sub(pbdot));
    var epsilon = 0.45;
    var num = -(1 + epsilon) * vrel;
    var term1 = 0;
    var term2 = 0;
    var term3 = 0;
    var term4 = 0;
    if (con.bodya.mass > 0)
    {
        term1 = 1 / con.bodya.mass;
        term3 = con.normal.dot(con.bodya.Iinv.vecMult(ra.cross(con.normal)).cross(ra));
    }
    if (con.bodyb.mass > 0)
    {
        term2 = 1 / con.bodyb.mass;
        term4 = con.normal.dot(con.bodyb.Iinv.vecMult(rb.cross(con.normal)).cross(rb));
    }
    var j = num / (term1 + term2 + term3 + term4);
    var impulse = con.normal.scalarMult(j);
    if (con.bodya.mass > 0)
    {
        con.bodya.P.accumulate(impulse);
        con.bodya.L.accumulate(ra.cross(impulse));
        con.bodya.applyAuxi();
    }
    if (con.bodyb.mass > 0)
    {
        con.bodyb.P.degeneratw(impulse);
        con.bodyb.L.degeneratw(rb.cross(impulse));
        con.bodyb.applyAuxi();
    }
}



PhysicsWorld.prototype.isColliding = function(con)
{
    var padot = this.ptVelocity(con.bodya, con.pos);
    var pbdot = this.ptVelocity(con.bodyb, con.pos);
    var vrel = con.normal.dot(padot.sub(pbdot));
    return (vrel < 0);
}

PhysicsWorld.prototype.isPenetrating = function(con)
{
    return (con.depth > 0.01);
}

PhysicsWorld.prototype.step = function (dt)
{
    //console.log("Step");
    
    var dt2 = dt / this.steps;
    for (var j = 0; j < this.steps; j++)
    {
        for (var i = 0; i < this.rigidBodies.length; i++)
        {
            this.rigidBodies[i].force.accumulate(this.gravity.scalarMult(this.rigidBodies[i].mass));
            // console.log("i: " + i+"-"+ this.rigidBodies[i].force.x+","+this.rigidBodies[i].force.y+","+this.rigidBodies[i].force.z);
        }
        this.copyBodies(this.prevRB, this.rigidBodies, false);
        this.integrateBodies(dt2);
        var contacts = this.findCollisionPoints();

        
        var collided;
        do
        {
            collided = false;
            for (var i = 0; i < contacts.length; i++)
            {
                var con = contacts[i];
                if (this.isColliding(con))
                {
                    this.handleCollision(con);
                }
            }
        }
        while (collided);
        for (var i = 0; i < this.rigidBodies.length; i++)
        {
            this.rigidBodies[i].force = new Vector(0, 0, 0);
            this.rigidBodies[i].torque = new Vector(0.1, 0.1, 0.1);
        }
    }

}
