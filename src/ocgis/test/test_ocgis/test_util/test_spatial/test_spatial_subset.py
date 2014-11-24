from copy import deepcopy

from shapely import wkt
import numpy as np

from ocgis import CoordinateReferenceSystem, RequestDataset
import ocgis
from ocgis.exc import EmptySubsetError
from ocgis.interface.base.crs import CFWGS84, CFRotatedPole, WrappableCoordinateReferenceSystem
from ocgis.interface.base.dimension.base import VectorDimension
from ocgis.interface.base.dimension.spatial import SpatialDimension
from ocgis.interface.base.field import Field
from ocgis.test.base import TestBase
from ocgis.test.test_ocgis.test_api.test_parms.test_definition import TestGeom
from ocgis.util.helpers import make_poly
from ocgis.util.itester import itr_products_keywords
from ocgis.util.spatial.spatial_subset import SpatialSubsetOperation
from ocgis import constants, env


class TestSpatialSubsetOperation(TestBase):

    def __init__(self, *args, **kwargs):
        self._target = None
        super(TestSpatialSubsetOperation, self).__init__(*args, **kwargs)

    def __iter__(self):
        keywords = dict(target=self.target,
                        output_crs=self.get_output_crs(),
                        wrap=[None, True, False])

        for k in itr_products_keywords(keywords, as_namedtuple=True):
            kwargs = k._asdict()
            target = kwargs.pop('target')
            ss = SpatialSubsetOperation(target, **kwargs)
            yield (ss, k)

    @property
    def germany(self):
        wkt_str = 'MULTIPOLYGON(((14.22597980499267578 53.92824935913085938,14.20333003997802557 53.90942001342773438,14.21741962432861328 53.86865997314453125,14.18721961975097656 53.87498092651367188,13.93943977355957031 53.84304046630859375,13.82499980926513672 53.85831832885742188,13.94054985046386541 53.91276931762695312,13.90305995941162109 53.98831939697265625,14.04889011383056463 53.94192886352539062,14.05249977111816406 53.99748992919921875,13.97663974761962891 54.04259109497070312,13.92000007629394531 54.06277084350585938,13.89054965972900391 54.00748825073241477,13.77554988861083984 54.02164840698242188,13.81721973419189453 54.10081863403319602,13.77554988861083984 54.13359832763671164,13.91139030456542969 54.08443069458007812,14.22597980499267578 53.92824935913085938)),((13.34082984924316406 54.23498916625976562,13.33555030822753906 54.27943038940428977,13.203610420227049 54.272491455078125,13.11833000183105469 54.33388137817382812,13.13082981109618963 54.370819091796875,13.26527976989746094 54.38026046752929688,13.15472030639648438 54.42443084716796164,13.27221965789794744 54.47581863403320312,13.14721965789794922 54.54082107543945312,13.24166011810302734 54.55276107788085938,13.30694007873534979 54.5133209228515625,13.36610984802246094 54.57526016235351562,13.34972000122070312 54.51942825317382102,13.37528038024902166 54.557769775390625,13.49860954284667791 54.47943115234375,13.51832962036132635 54.56248092651367188,13.40167045593261541 54.57304000854492188,13.36528015136718572 54.61526107788085938,13.28639030456542969 54.56332015991210938,13.24610996246337891 54.557769775390625,13.28888988494872869 54.64471054077147727,13.22665977478027344 54.62860107421875,13.23054981231689275 54.64776992797851562,13.41444015502929688 54.6819305419921875,13.38304996490478516 54.63888168334960938,13.4002799987792951 54.59415054321289062,13.45499992370605469 54.57276916503905539,13.63389015197753906 54.58581924438475852,13.67667007446289062 54.56526947021483664,13.67527961730957031 54.52637100219726562,13.58528041839599609 54.48276901245117188,13.57777976989746094 54.45388031005859375,13.61750030517578125 54.40414810180664062,13.76778030395507635 54.3408203125,13.69639015197753729 54.29275894165038352,13.69054985046386719 54.32720947265625,13.68443965911865234 54.349151611328125,13.44972038269042791 54.31637954711914062,13.36499977111816406 54.26137924194335938,13.42000007629394531 54.26137924194335938,13.42943954467773438 54.23804092407225852,13.34082984924316406 54.23498916625976562)),((9.44516944885253906 54.82477951049804688,9.50312423706054688 54.8482513427734375,9.45166587829589844 54.8069305419921875,9.57694435119628906 54.86193084716796875,9.57833290100097656 54.82555007934569602,9.91666603088378906 54.79109954833983664,10.01583003997802734 54.6955413818359375,9.93093109130859375 54.66756057739257102,10.03777980804443359 54.6683197021484375,10.02777957916259766 54.55110168457030539,9.86583328247070312 54.47748947143554688,10.12693977355957031 54.48971176147460938,10.20277976989745916 54.45804977416992188,10.19806003570556641 54.39582061767578125,10.15167045593261719 54.3636016845703125,10.23805046081542969 54.41720962524414062,10.37250041961669922 54.4341583251953125,10.66499996185302557 54.32442855834960938,10.77999973297119141 54.307769775390625,10.84889030456542969 54.34193038940429688,11.12222003936767578 54.39249038696288352,11.12917041778564453 54.37554168701171875,11.06472015380859375 54.3499908447265625,11.09305000305175781 54.19749069213867188,10.87388992309570312 54.08776855468749289,10.80222034454345703 54.09553909301757812,10.75638961791992188 54.03443145751953125,10.78999996185302734 53.99721145629882812,10.87055015563964844 53.99193954467773438,10.90262985229492188 53.95998001098632102,11.04582977294921875 54.00664901733397727,11.19221973419189453 54.01054000854492188,11.24221992492675781 53.94470977783203125,11.27694034576416016 53.93165969848632812,11.33804988861083807 53.95721054077148438,11.45417022705078125 53.90110015869140625,11.4713897705078125 53.96776962280273438,11.58638954162597479 54.04248809814453125,11.62444019317626953 54.11027145385742188,12.02000045776367188 54.17971038818358664,12.10000038146972656 54.18165969848632812,12.34360980987548828 54.29777145385741477,12.48721981048583984 54.45499038696289062,12.91166019439697266 54.43970870971678977,12.91749954223632812 54.41999053955078125,12.77583026885986328 54.41666030883789062,12.60583019256591797 54.40610122680664062,12.46000003814697266 54.39220809936523438,12.369720458984375 54.26499176025389914,12.41804981231689453 54.2608184814453125,12.42055988311767578 54.28525924682616477,12.46193981170654297 54.302490234375,12.71249961853027344 54.40943145751953125,12.74722003936767578 54.37916183471679688,12.81610965728759766 54.34637069702148438,12.99806022644042791 54.43442916870117188,13.0922203063964826 54.3694305419921875,13.07888984680175781 54.33166122436522727,13.16082954406738104 54.26248931884765625,13.38833045959472479 54.14720916748046164,13.49971961975097656 54.08610153198242188,13.711669921875 54.17137908935546875,13.80805015563964844 54.09693145751953125,13.75055980682373047 54.02859878540039062,13.91306018829345703 53.91748809814453125,13.80749988555908203 53.85525894165039062,13.9036102294921875 53.80276107788085938,14.14083003997802734 53.7386016845703125,14.23639011383056641 53.75915908813476562,14.25771045684814453 53.7431793212890625,14.23750019073486328 53.69860076904296875,14.27781963348388672 53.69392013549804688,14.28639030456542791 53.66915130615234375,14.37555027008056641 53.42303848266601562,14.44530010223388494 53.27259063720702414,14.40861034393310547 53.2158203125,14.39083003997802557 53.14165115356445312,14.3422203063964826 53.04499053955078125,14.16471958160400391 52.96886825561522727,14.13333034515380682 52.83332061767578125,14.25749969482421875 52.79053878784179688,14.55111026763916016 52.62858963012695312,14.64111042022705078 52.56666183471679688,14.53444004058837891 52.39471054077148438,14.5669403076171875 52.32804107666015625,14.7133302688598615 52.2391510009765625,14.69305038452148438 52.10470962524413352,14.76352977752685369 52.07080841064452414,14.747220039367674 52.05636978149414062,14.71638965606689631 51.94108963012695312,14.61610984802245916 51.85387039184570312,14.6019401550292951 51.8135986328125,14.75749969482421875 51.65942001342773438,14.71749973297118963 51.55276107788085938,14.91082954406738459 51.48303985595703125,15.03804969787597479 51.26803970336914062,14.97527980804443359 51.10638046264647727,14.82884979248046875 50.86603164672851562,14.80165958404541016 50.818878173828125,14.62526988983154297 50.85416030883788352,14.65083026885986328 50.92416000366210938,14.56972026824950994 50.91608810424804688,14.59889030456542791 50.97719955444335938,14.47222042083740234 51.03136825561523438,14.26082992553710938 50.99665069580078125,14.4036102294921875 50.93249130249022727,14.39361000061035156 50.89471054077148438,14.04666996002197266 50.80693054199218039,13.90777969360351385 50.79053878784179688,13.86194038391113281 50.72581100463866477,13.56667041778564275 50.71110153198242188,13.47000026702880682 50.60332107543945312,13.38000011444091797 50.64136886596678977,13.32083034515380859 50.58110046386718039,13.23943996429443359 50.57970809936523438,13.19110965728759588 50.50360107421875,13.03361034393310547 50.50276947021484375,12.99444007873534979 50.42303848266601562,12.90388965606689453 50.42332077026367188,12.7849998474121076 50.44636917114257812,12.48639011383056641 50.35137939453125,12.36472034454345703 50.27914810180664062,12.32971954345702947 50.16970825195311789,12.28055000305175781 50.18610000610351562,12.26000022888183594 50.26166152954101562,12.09346961975097479 50.32419967651366477,12.13582992553710938 50.2783203125,12.09694004058837891 50.2497100830078125,12.20499992370605469 50.17416000366210938,12.22889041900634766 50.09637069702148438,12.42889022827148438 49.98442840576171875,12.49277019500732422 49.97497940063475852,12.47443962097167969 49.94303894042968039,12.54249954223632635 49.92026901245117188,12.49944019317626953 49.83275985717773438,12.40305042266845703 49.75999069213867188,12.44194030761718572 49.70026016235351562,12.52499961853027344 49.63721084594726562,12.66166019439697088 49.433868408203125,12.87860965728759766 49.32804107666015625,12.93527984619140625 49.34025955200195312,13.02694034576416016 49.29832077026367188,13.15110969543456854 49.17747879028319602,13.203610420227049 49.1194305419921875,13.39527988433837891 49.05025863647460938,13.4014396667480451 49.00748825073242188,13.50166988372802734 48.94498062133789062,13.51611042022705078 48.97776031494139914,13.58444023132324041 48.96886825561523438,13.66193962097167969 48.89638137817382812,13.74083042144775213 48.88164901733398438,13.81462955474853516 48.787139892578125,13.79304981231689453 48.7260894775390625,13.83333015441894531 48.69887161254882812,13.72665977478027344 48.51776123046874289,13.50527954101562322 48.58304977416992188,13.4483299255371076 48.56859970092773438,13.43389034271240234 48.41997909545898438,13.36861038208007812 48.3519287109375,12.9295196533203125 48.20930099487304688,12.75666999816894531 48.12054061889648438,13.00889015197753906 47.85416030883788352,12.93972015380859375 47.78470993041992188,12.91749954223632812 47.71554183959960938,13.05832958221435547 47.70608901977539062,13.10027980804443359 47.64081954956053977,13.07417011260986328 47.6169281005859375,13.05305004119873047 47.49637985229492188,13.00500011444091797 47.46942901611328125,12.809720039367674 47.55221176147460938,12.78777980804443359 47.58942031860351562,12.83082962036132812 47.61886978149414062,12.77388954162597479 47.67416000366210227,12.50500011444091797 47.63748931884765625,12.44167041778564453 47.69858932495117188,12.24388980865478516 47.69469833374022727,12.25722026824951172 47.74303817749023438,12.17722034454345703 47.70164871215819602,12.206390380859375 47.64136886596678977,12.19610977172851562 47.60942840576171875,11.87860965728759766 47.60665130615233664,11.63333034515380859 47.5952606201171875,11.57444000244140625 47.51998138427733664,11.42193984985351562 47.508880615234375,11.40555000305175781 47.45386886596679688,11.22749996185302734 47.40053939819335938,11.2366600036621076 47.43304061889648438,11.20304965972900391 47.43526077270507812,11.10694026947021307 47.39638137817382812,10.97360992431640625 47.40053939819335938,10.95028018951416016 47.46026992797851562,10.86583042144775391 47.49303817749023438,10.90972042083740234 47.52193069458007812,10.84860992431640625 47.53638076782226562,10.68554973602294922 47.55858993530273438,10.55443954467773438 47.53693008422851562,10.46527004241943359 47.55858993530273438,10.48276996612548828 47.59054183959960938,10.42693996429443182 47.57693099975585938,10.47332954406738281 47.43553924560546164,10.27443981170654297 47.28887939453124289,10.16971969604492188 47.28110122680663352,10.21415996551513672 47.3152618408203125,10.21138954162597656 47.38637924194335938,10.15443992614746094 47.36914825439452414,10.08722019195556641 47.38721084594725852,10.10443973541259588 47.42887115478515625,10.08444023132324219 47.46026992797851562,10.00333023071289062 47.4838714599609375,9.96333122253417969 47.54777145385742188,9.85416412353515625 47.53887939453125,9.81361007690429688 47.59360122680663352,9.76305389404296875 47.58470916748046875,9.72729301452636719 47.53625869750975852,9.60555267333984375 47.52914810180663352,9.5676116943359375 47.54391860961913352,9.26166534423828125 47.6630401611328125,8.87860870361328125 47.65581893920898438,8.79999923706054688 47.735260009765625,8.79749870300292969 47.68304061889648438,8.72499847412109375 47.69776153564453125,8.72694206237792969 47.76499176025390625,8.55972099304199219 47.80636978149413352,8.40666580200195312 47.70386886596679688,8.41333198547363281 47.67110061645507812,8.62166404724120916 47.66025924682617188,8.61888694763183594 47.63970947265625,8.50833320617675781 47.62831878662108664,8.49303054809570312 47.58456039428710938,8.34055328369140625 47.57416152954100852,8.20555305480957031 47.62165069580078125,7.94333219528198242 47.55360031127929688,7.81944417953491122 47.58832168579101562,7.69722080230712891 47.5433197021484375,7.6183319091796875 47.56110000610350852,7.67444276809692294 47.60638046264648438,7.58879899978637695 47.58456039428710938,7.51333284378051758 47.68693161010741477,7.61583185195922852 48.00276947021484375,7.57166576385497958 48.03720855712890625,7.60499906539916903 48.15692901611328125,7.75083303451537997 48.33665084838867188,7.7716660499572745 48.49164962768553977,7.80749893188476474 48.5133209228515625,7.80194282531738281 48.59247970581053977,7.92151498794555575 48.69002914428710938,8.13333320617675781 48.88554000854492188,8.22739410400390625 48.96371078491210938,8.19222068786621094 48.96886825561523438,7.93861007690429688 49.04887008666992188,7.74110984802246094 49.04166030883789062,7.53999900817871094 49.08887100219726562,7.48694276809692383 49.16415023803710227,7.37527704238891513 49.17192840576171164,7.36138820648193359 49.14776992797851562,7.08833217620849609 49.1252593994140625,7.03833198547363281 49.11832046508789062,7.02666616439819336 49.18886947631835938,6.84749889373779297 49.21525955200195312,6.83888816833496094 49.15497970581054688,6.78722000122070312 49.16247940063476562,6.58972215652465731 49.32027053833007812,6.59472179412841797 49.36304092407225852,6.49388790130615234 49.44720077514648438,6.36222076416015625 49.45998001098632812,6.37249898910522461 49.59025955200194602,6.51055479049682617 49.70637893676757812,6.52222204208374023 49.81110000610351562,6.32666587829589844 49.83971023559570312,6.13972187042236328 49.99665069580078125,6.13183307647705078 50.12553024291992188,6.17333221435546875 50.23247909545898438,6.40027713775634677 50.32915115356445312,6.36638689041137695 50.45220947265625,6.20083189010620028 50.51638031005859375,6.26861000061035067 50.62360000610351562,6.17138814926147461 50.623870849609375,6.10833311080932617 50.72330856323242188,6.0286102294921875 50.7158203125,6.0084071159362793 50.75606918334960938,5.98166608810424805 50.80276107788085938,6.08472204208374023 50.87360000610350852,6.01083278656005859 50.94359970092773438,6.02833318710327148 50.97665023803710938,5.96583318710327148 50.97832107543945312,5.86944293975830078 51.01887893676757812,5.87388801574707031 51.05025863647460938,5.95222187042236239 51.036651611328125,6.16722106933593661 51.16276168823242188,6.07944393157958984 51.17581939697265625,6.07388877868652344 51.22053909301757812,6.22222089767455966 51.36166000366210938,6.217498779296875 51.47637176513671875,6.09305477142333984 51.607208251953125,6.11610984802246094 51.65192031860351562,6.02944421768188477 51.67860031127928977,6.03916501998901367 51.71693038940429688,5.95499897003173739 51.73859024047851562,5.9688878059387207 51.79109954833983664,5.9619441032409668 51.83026885986328125,6.1697211265563956 51.84193038940429688,6.13777685165405273 51.87693023681640625,6.15972089767456055 51.90554046630859375,6.3802771568298331 51.82999038696289062,6.54888820648193359 51.88526153564453125,6.72777700424194336 51.8994293212890625,6.83083200454711914 51.97137069702148438,6.80055379867553622 52.00720977783203125,6.68805408477783203 52.03887939453125,6.69749784469604492 52.06998062133788352,6.86055517196655185 52.12025833129882812,6.87972116470336825 52.15359115600585938,7.05249881744384677 52.23582077026367188,7.06555509567260742 52.3858184814453125,6.98749923706054599 52.46110153198242188,6.94610977172851562 52.43442916870117188,6.70555496215820312 52.48582077026367188,6.69027614593505859 52.55192947387695312,6.76083278656005859 52.56721115112304688,6.72083282470703125 52.62942886352538352,6.78166580200195312 52.65414810180664062,7.03583192825317383 52.63275909423828125,7.05527687072753906 52.65192031860351562,7.06972122192382812 52.81499099731444602,7.19888782501220703 52.9677581787109375,7.20944404602050781 53.24275970458984375,7.21111106872558594 53.24415969848632812,7.25277709960937589 53.31692886352539062,7.05055522918701172 53.33971023559570312,7.01583290100097567 53.38359832763671875,7.03499984741210938 53.48749160766601562,7.14722204208374023 53.53720855712890625,7.08888816833496005 53.57109832763671875,7.09610986709594727 53.59193038940429688,7.25250005722045898 53.67387008666991477,7.39083290100097567 53.68637847900390625,7.84861087799072266 53.71416091918944602,8.02638816833496094 53.70360183715820312,8.05055427551269531 53.63220977783203125,8.15944290161132812 53.56248092651367188,8.16333198547363281 53.52054977416992188,8.06138801574707031 53.50054168701171875,8.08027553558349432 53.45804977416992188,8.22305488586425781 53.40082931518554688,8.28416633605957031 53.41859817504882812,8.3138885498046875 53.45943069458007812,8.31194305419921875 53.52415084838867188,8.23194313049316406 53.52220916748046875,8.24305534362792969 53.58610153198241477,8.28333282470702947 53.61388015747069602,8.34138870239257812 53.61360168457030539,8.54888725280761719 53.52943038940429688,8.48833274841308594 53.47916030883789062,8.48527717590332031 53.40610122680664062,8.502777099609375 53.40887069702148438,8.56833267211914062 53.53332901000976562,8.51972198486328125 53.60303878784179688,8.4877777099609375 53.70193099975585938,8.60416603088378906 53.87942886352539062,8.68055534362792969 53.89471054077148438,8.72805404663085938 53.85776901245117188,8.90749931335449219 53.82804107666015625,9.12888717651367188 53.86610031127928977,9.29521560668945312 53.83568954467773438,9.35999870300292969 53.78638076782225852,9.58277702331542969 53.58610153198241477,9.80222129821777344 53.53470993041991477,9.82472038269042791 53.55027008056640625,9.67416572570800781 53.57442855834960227,9.54666519165039062 53.63359832763671875,9.52972221374511719 53.70775985717773438,9.43777656555175781 53.73942947387695312,9.3763885498046875 53.83137893676757812,8.96194076538085938 53.89860153198241477,8.85694313049316406 54.00276947021484375,8.85777664184570135 54.04220962524414062,9.00305366516113281 54.02804946899414062,9.01777553558349609 54.09109878540039062,8.97027587890625 54.14554977416992188,8.84472084045410156 54.1330413818359375,8.82333183288574219 54.20721054077148438,8.91611099243164062 54.27082061767578125,8.95805454254150391 54.31692886352539062,8.79388809204101562 54.28416061401367188,8.63166618347167969 54.27914810180664062,8.5994415283203125 54.32749176025390625,8.60916519165039062 54.34553909301757102,8.69388771057128906 54.35749053955078125,8.65416526794433416 54.37582015991210938,8.90166664123535156 54.41942977905272727,9.01472091674804688 54.48025894165039062,8.97495841979980469 54.51433181762695312,8.85083198547363104 54.620819091796875,8.80562210083007812 54.68521881103515625,8.67610931396484375 54.77914810180663352,8.64222145080566406 54.82638168334960938,8.6562957763671875 54.91740036010742188,8.71944236755371094 54.89110183715820312,8.91999816894531072 54.90803909301757812,9.44516944885253906 54.82477951049804688)))'
        germany = self.get_buffered(wkt.loads(wkt_str))
        germany = {'geom': germany, 'properties': {'UGID': 2, 'DESC': 'Germany'}}
        return germany

    @property
    def nebraska(self):
        wkt_str = 'POLYGON((-101.407393 40.001003,-102.051535 39.998918,-102.047545 40.342644,-102.047620 40.431077,-102.046031 40.697319,-102.046992 40.743130,-102.047739 40.998071,-102.621257 41.000214,-102.652271 40.998124,-103.382956 41.000316,-103.572316 40.999648,-104.051705 41.003211,-104.054012 41.388085,-104.055500 41.564222,-104.053615 41.698218,-104.053513 41.999815,-104.056219 42.614669,-104.056199 43.003062,-103.501464 42.998618,-103.005875 42.999354,-102.788384 42.995303,-102.086701 42.989887,-101.231737 42.986843,-100.198142 42.991095,-99.532790 42.992335,-99.253971 42.992389,-98.497651 42.991778,-98.457444 42.937160,-98.391204 42.920135,-98.310339 42.881794,-98.167826 42.839571,-98.144869 42.835794,-98.123117 42.820223,-98.121820 42.808360,-98.033140 42.769192,-97.995144 42.766812,-97.963558 42.773690,-97.929477 42.792324,-97.889941 42.831271,-97.888659 42.855807,-97.818643 42.866587,-97.797028 42.849597,-97.772186 42.846164,-97.725250 42.858008,-97.685752 42.836837,-97.634970 42.861285,-97.570654 42.847990,-97.506132 42.860136,-97.483159 42.857157,-97.457263 42.850443,-97.389306 42.867433,-97.311414 42.861771,-97.271457 42.850014,-97.243189 42.851826,-97.224443 42.841202,-97.211831 42.812573,-97.161422 42.798619,-97.130469 42.773923,-97.015139 42.759542,-96.979593 42.758313,-96.970003 42.752065,-96.977869 42.727308,-96.970773 42.721147,-96.908234 42.731699,-96.810140 42.704084,-96.810437 42.681341,-96.799344 42.670019,-96.722658 42.668592,-96.699060 42.657715,-96.694596 42.641163,-96.715273 42.621907,-96.714059 42.612302,-96.636672 42.550731,-96.629294 42.522693,-96.605467 42.507236,-96.584753 42.518287,-96.547215 42.520499,-96.494701 42.488459,-96.439394 42.489240,-96.396074 42.467401,-96.397890 42.441793,-96.417628 42.414777,-96.411761 42.380918,-96.424175 42.349279,-96.389781 42.328789,-96.368700 42.298023,-96.342881 42.282081,-96.332658 42.260307,-96.337708 42.229522,-96.363512 42.214042,-96.352165 42.168185,-96.285123 42.123452,-96.265483 42.048897,-96.238725 42.028438,-96.236093 42.001258,-96.202842 41.996615,-96.185217 41.980685,-96.147328 41.966254,-96.145870 41.924907,-96.159970 41.904151,-96.135623 41.862620,-96.076417 41.791469,-96.099321 41.752975,-96.099771 41.731563,-96.085557 41.704987,-96.122202 41.694913,-96.120264 41.684094,-96.099306 41.654680,-96.111307 41.599006,-96.080835 41.576000,-96.091936 41.563145,-96.085840 41.537522,-96.050172 41.524335,-96.004592 41.536663,-95.993965 41.528103,-95.996688 41.511517,-96.013451 41.492994,-96.006897 41.481954,-95.953185 41.472387,-95.935065 41.462381,-95.940056 41.394805,-95.942895 41.340077,-95.889107 41.301389,-95.897591 41.286863,-95.911202 41.308469,-95.930230 41.302056,-95.910981 41.225245,-95.922250 41.207854,-95.916100 41.194063,-95.859198 41.180537,-95.859801 41.166865,-95.876685 41.164202,-95.858274 41.109187,-95.878804 41.065871,-95.859539 41.035002,-95.860897 41.002650,-95.837603 40.974258,-95.836541 40.901108,-95.834396 40.870300,-95.846435 40.848332,-95.851790 40.792600,-95.876616 40.730436,-95.767999 40.643117,-95.757546 40.620904,-95.767479 40.589048,-95.763412 40.549707,-95.737036 40.532373,-95.692066 40.524129,-95.687413 40.561170,-95.675693 40.565835,-95.662944 40.558729,-95.658060 40.530332,-95.684970 40.512205,-95.695361 40.485338,-95.636817 40.396390,-95.634185 40.358800,-95.616201 40.346497,-95.617933 40.331418,-95.645553 40.322346,-95.646827 40.309109,-95.595532 40.309776,-95.547137 40.266215,-95.476822 40.226855,-95.466636 40.213255,-95.460952 40.173995,-95.422476 40.131743,-95.392813 40.115416,-95.384542 40.095362,-95.403784 40.080379,-95.413764 40.048111,-95.390532 40.043750,-95.371244 40.028751,-95.345067 40.024974,-95.308697 39.999407,-95.329701 39.992595,-95.780700 39.993489,-96.001253 39.995159,-96.240598 39.994503,-96.454038 39.994172,-96.801420 39.994476,-96.908287 39.996154,-97.361912 39.997380,-97.816589 39.999729,-97.929588 39.998452,-98.264165 39.998434,-98.504479 39.997129,-98.720632 39.998461,-99.064747 39.998338,-99.178201 39.999577,-99.627859 40.002987,-100.180910 40.000478,-100.191111 40.000585,-100.735049 39.999172,-100.754856 40.000198,-101.322148 40.001821,-101.407393 40.001003))'
        nebraska = self.get_buffered(wkt.loads(wkt_str))
        nebraska = {'geom': nebraska, 'properties': {'UGID': 1, 'DESC': 'Nebraska'}}
        return nebraska

    @property
    def rd_rotated_pole(self):
        rd = self.test_data.get_rd('rotated_pole_cccma')
        return rd

    @property
    def target(self):
        if self._target is None:
            self._target = self.get_target()
        return self._target

    def get_buffered(self, geom):
        ret = geom.buffer(0)
        self.assertTrue(ret.is_valid)
        return ret

    def get_output_crs(self):
        crs_wgs84 = CFWGS84()
        ret = ['input', crs_wgs84]
        return ret

    def get_subset_sdim(self):

        # 1: nebraska
        nebraska = self.nebraska

        # 2: germany
        germany = self.germany

        # 3: nebraska and germany

        ret = [SpatialDimension.from_records(d) for d in [[nebraska], [germany], [nebraska, germany]]]

        return ret

    def get_target(self):
        # 1: standard input file - geographic coordinate system, unwrapped
        rd_standard = self.test_data.get_rd('cancm4_tas')

        # 2: standard field - geographic coordinate system
        field_standard = rd_standard.get()

        # 3: field with rotated pole coordinate system
        field_rotated_pole = self.rd_rotated_pole.get()

        # 4: field with lambert conformal coordinate system
        rd = self.test_data.get_rd('narccap_lambert_conformal')
        field_lambert = rd.get()

        # 5: standard input field - geographic coordinate system, wrapped
        field_wrapped = rd_standard.get()
        field_wrapped.spatial.wrap()

        # 6: spatial dimension - standard geographic coordinate system
        sdim = field_standard.spatial

        ret = [rd_standard, field_standard, field_rotated_pole, field_lambert, field_wrapped, sdim]

        return ret

    def test_init_output_crs(self):
        for ss, k in self:
            if k.output_crs is None:
                if isinstance(k.target, Field):
                    self.assertEqual(ss.sdim.crs, k.target.spatial.crs)

    def test_field(self):
        for ss, k in self:
            try:
                self.assertIsInstance(ss.field, Field)
            except AttributeError:
                if isinstance(k.target, SpatialDimension):
                    continue
                else:
                    raise

    def test_get_buffered_subset_sdim(self):
        proj4 = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'
        buffer_crs_list = [None, CoordinateReferenceSystem(proj4=proj4)]
        poly = make_poly((36, 44), (-104, -95))

        for buffer_crs in buffer_crs_list:
            subset_sdim = SpatialDimension.from_records([{'geom': poly, 'properties': {'UGID': 1}}], crs=CFWGS84())
            self.assertEqual(subset_sdim.crs, CFWGS84())
            if buffer_crs is None:
                buffer_value = 1
            else:
                buffer_value = 10

            ret = SpatialSubsetOperation._get_buffered_subset_sdim_(subset_sdim, buffer_value, buffer_crs=buffer_crs)
            ref = ret.geom.polygon.value[0, 0]

            if buffer_crs is None:
                self.assertEqual(ref.bounds, (-105.0, 35.0, -94.0, 45.0))
            else:
                self.assertNumpyAllClose(np.array(ref.bounds), np.array((-104.00013263459613, 35.9999147913708, -94.99986736540386, 44.00008450528758)))
            self.assertEqual(subset_sdim.crs, ret.crs)

            # check deepcopy
            ret.geom.polygon.value[0, 0] = make_poly((1, 2), (3, 4))
            ref_buffered = ret.geom.polygon.value[0, 0]
            ref_original = subset_sdim.geom.polygon.value[0, 0]
            with self.assertRaises(AssertionError):
                self.assertNumpyAllClose(np.array(ref_buffered.bounds), np.array(ref_original.bounds))

    def test_get_should_wrap(self):
        # a 360 dataset
        field_360 = self.test_data.get_rd('cancm4_tas').get()
        ss = SpatialSubsetOperation(field_360, wrap=True)
        self.assertTrue(ss._get_should_wrap_(ss.target))
        ss = SpatialSubsetOperation(field_360, wrap=False)
        self.assertFalse(ss._get_should_wrap_(ss.target))
        ss = SpatialSubsetOperation(field_360, wrap=None)
        self.assertFalse(ss._get_should_wrap_(ss.target))

        # wrapped dataset
        field_360.spatial.wrap()
        ss = SpatialSubsetOperation(field_360, wrap=True)
        self.assertFalse(ss._get_should_wrap_(ss.target))
        ss = SpatialSubsetOperation(field_360, wrap=False)
        self.assertFalse(ss._get_should_wrap_(ss.target))
        ss = SpatialSubsetOperation(field_360, wrap=None)
        self.assertFalse(ss._get_should_wrap_(ss.target))

    def test_get_spatial_subset(self):
        ctr_test = 0
        ctr = 0
        for ss, k in self:
            for subset_sdim in self.get_subset_sdim():
                for operation in ['intersects', 'clip', 'foo']:

                    use_subset_sdim = deepcopy(subset_sdim)
                    use_ss = deepcopy(ss)

                    # ctr += 1
                    # print ctr
                    # if ctr != 73:
                    #     continue
                    # else:
                    #     import ipdb;ipdb.set_trace()

                    try:
                        ret = use_ss.get_spatial_subset(operation, use_subset_sdim, use_spatial_index=True,
                                                        select_nearest=False, buffer_value=None, buffer_crs=None)
                    except ValueError:
                        # 'foo' is not a valid type of subset operation.
                        if operation == 'foo':
                            continue
                        # only one polygon for a spatial operation
                        elif use_subset_sdim.shape != (1, 1):
                            continue
                        else:
                            raise
                    except EmptySubsetError:
                        # subset tests occur on the spatial dimension operations
                        continue
                    try:
                        self.assertIsInstance(ret, type(use_ss.target))
                    except AssertionError:
                        # if the target is a request datasets, then the output should be a field
                        if isinstance(use_ss.target, RequestDataset):
                            self.assertIsInstance(ret, Field)
                    ctr_test += 1
        self.assertGreater(ctr_test, 5)

    def test_get_spatial_subset_circular_geometries(self):
        """Test circular geometries. They were causing wrapping errors."""

        geoms = TestGeom.get_geometry_dictionaries()
        rd = self.test_data.get_rd('cancm4_tas')
        ss = SpatialSubsetOperation(rd, wrap=True)
        buffered = [element['geom'].buffer(rd.get().spatial.grid.resolution*2) for element in geoms]
        for buff in buffered:
            record = [{'geom': buff, 'properties': {'UGID': 1}}]
            subset_sdim = SpatialDimension.from_records(record)
            ret = ss.get_spatial_subset('intersects', subset_sdim)
            self.assertTrue(np.all(ret.spatial.grid.extent > 0))

    def test_get_spatial_subset_output_crs(self):
        """Test subsetting with an output CRS."""

        # test with default crs converting to north american lambert
        proj4 = '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'
        output_crs = CoordinateReferenceSystem(proj4=proj4)
        subset_sdim = SpatialDimension.from_records([self.nebraska])
        rd = self.test_data.get_rd('cancm4_tas')
        ss = SpatialSubsetOperation(rd, output_crs=output_crs)
        ret = ss.get_spatial_subset('intersects', subset_sdim)
        self.assertEqual(ret.spatial.crs, output_crs)
        self.assertAlmostEqual(ret.spatial.grid.value.mean(), -35065.750850951554)

        # test with an input rotated pole coordinate system
        rd = self.rd_rotated_pole
        ss = SpatialSubsetOperation(rd, output_crs=env.DEFAULT_COORDSYS)
        subset_sdim = SpatialDimension.from_records([self.germany])
        ret = ss.get_spatial_subset('intersects', subset_sdim)
        self.assertEqual(ret.spatial.crs, env.DEFAULT_COORDSYS)

    def test_get_spatial_subset_rotated_pole(self):
        """Test input has rotated pole with now output CRS."""

        rd = self.rd_rotated_pole
        ss = SpatialSubsetOperation(rd)
        subset_sdim = SpatialDimension.from_records([self.germany])
        ret = ss.get_spatial_subset('intersects', subset_sdim)
        self.assertEqual(ret.spatial.crs, rd.get().spatial.crs)
        self.assertAlmostEqual(ret.spatial.grid.value.data.mean(), -2.0600000000000009)

    def test_get_spatial_subset_wrap(self):
        """Test subsetting with wrap set to a boolean value."""

        subset_sdim = SpatialDimension.from_records([self.nebraska])
        rd = self.test_data.get_rd('cancm4_tas')
        self.assertEqual(rd.get().spatial.wrapped_state, WrappableCoordinateReferenceSystem._flag_unwrapped)
        ss = SpatialSubsetOperation(rd, wrap=True)
        ret = ss.get_spatial_subset('intersects', subset_sdim)
        self.assertEqual(ret.spatial.wrapped_state, WrappableCoordinateReferenceSystem._flag_wrapped)
        self.assertAlmostEqual(ret.spatial.grid.value.data[1].mean(), -99.84375)

        # test with wrap false
        ss = SpatialSubsetOperation(rd, wrap=False)
        ret = ss.get_spatial_subset('intersects', subset_sdim)
        self.assertEqual(ret.spatial.wrapped_state, WrappableCoordinateReferenceSystem._flag_unwrapped)
        self.assertAlmostEqual(ret.spatial.grid.value.data[1].mean(), 260.15625)

    def test_prepare_target(self):
        for ss, k in self:
            self.assertIsNone(ss._original_rotated_pole_state)
            if isinstance(ss.sdim.crs, CFRotatedPole):
                ss._prepare_target_()
                self.assertIsInstance(ss._original_rotated_pole_state, CFRotatedPole)
                self.assertIsInstance(ss.sdim.crs, CFWGS84)
            else:
                ss._prepare_target_()
                self.assertIsNone(ss._original_rotated_pole_state)

    def test_prepare_subset_sdim(self):
        for subset_sdim in self.get_subset_sdim():
            for ss, k in self:
                try:
                    prepared = ss._prepare_subset_sdim_(subset_sdim)
                    # check that a deepcopy has occurred
                    self.assertFalse(np.may_share_memory(prepared.uid, subset_sdim.uid))
                except KeyError:
                    # the target has a rotated pole coordinate system. transformations to rotated pole for the subset
                    # geometry is not supported.
                    if isinstance(ss.sdim.crs, CFRotatedPole):
                        continue
                    else:
                        raise
                self.assertEqual(prepared.crs, ss.sdim.crs)

        # test nebraska against an unwrapped dataset specifically
        nebraska = SpatialDimension.from_records([self.nebraska])
        field = self.test_data.get_rd('cancm4_tas').get()
        ss = SpatialSubsetOperation(field)
        prepared = ss._prepare_subset_sdim_(nebraska)
        self.assertEqual(prepared.wrapped_state, WrappableCoordinateReferenceSystem._flag_unwrapped)

    def test_sdim(self):
        for ss, k in self:
            self.assertIsInstance(ss.sdim, SpatialDimension)

    def test_should_update_crs(self):
        # no output crs provided
        target = self.test_data.get_rd('cancm4_tas')
        ss = SpatialSubsetOperation(target)
        self.assertFalse(ss.should_update_crs)

        # output crs different than input
        ss = SpatialSubsetOperation(target, output_crs=CoordinateReferenceSystem(epsg=2136))
        self.assertTrue(ss.should_update_crs)

        # same output crs as input
        ss = SpatialSubsetOperation(target, output_crs=ss.sdim.crs)
        self.assertFalse(ss.should_update_crs)
