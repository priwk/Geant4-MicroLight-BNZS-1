#ifndef StageDOpticalEventAction_h
#define StageDOpticalEventAction_h 1

#include "G4UserEventAction.hh"
#include "StageDOpticalStats.hh"

class G4Event;
class AnalysisConfig;
class StageDOpticalPrimaryGeneratorAction;
class StageDOpticalRunAction;

class StageDOpticalEventAction : public G4UserEventAction
{
public:
  StageDOpticalEventAction(StageDOpticalRunAction *runAction,
                           const StageDOpticalPrimaryGeneratorAction *primaryAction,
                           AnalysisConfig *config);
  ~StageDOpticalEventAction() override;

  void BeginOfEventAction(const G4Event *event) override;
  void EndOfEventAction(const G4Event *event) override;

  StageDPhotonEventRecord &MutableCurrentEvent() { return fCurrentEvent; }
  const StageDPhotonEventRecord &GetCurrentEvent() const { return fCurrentEvent; }

  void SetFinalStatus(const std::string &status, G4bool absorbed);
  void MarkAbsorbed(const std::string &phaseLabel);
  void MarkCensoredEncounterIfActive(const std::string &reason);

private:
  AnalysisConfig *fConfig;
  StageDOpticalRunAction *fRunAction;
  const StageDOpticalPrimaryGeneratorAction *fPrimaryAction;
  StageDPhotonEventRecord fCurrentEvent;
};

#endif
